use std::path::Path;

use rayon::prelude::*;
use rsomics_common::Result;
use rsomics_common::RsomicsError;
use rsomics_fqgz::ChunkedWriter;
use rsomics_seqio::{OwnedRecord, open_fastq};
use serde::Serialize;

use crate::extract::process;
use crate::umi_loc::UmiConfig;

const CHUNK_RECORDS: usize = 8192;

#[derive(Debug, Default, Clone, Serialize)]
pub struct UmiReport {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mode: Option<&'static str>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub input_r1: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub input_r2: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub output_r1: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub output_r2: Option<String>,
    pub reads_in: u64,
    pub reads_out: u64,
    pub bases_in: u64,
    pub bases_out: u64,
}

struct OwnedPair {
    r1: OwnedRecord,
    r2: OwnedRecord,
}

pub struct Pipeline<'cfg> {
    pub cfg: &'cfg UmiConfig,
    pub compression: i32,
}

impl<'cfg> Pipeline<'cfg> {
    #[must_use]
    pub fn new(cfg: &'cfg UmiConfig, compression: i32) -> Self {
        Self { cfg, compression }
    }

    pub fn run_se(&self, input: &Path, output: &Path) -> Result<UmiReport> {
        let mut reader = open_fastq(input)?;
        let mut writer = ChunkedWriter::create(output, self.compression)?;
        let mut report = UmiReport {
            mode: Some("SE"),
            input_r1: Some(input.display().to_string()),
            output_r1: Some(output.display().to_string()),
            ..UmiReport::default()
        };
        let mut chunk: Vec<OwnedRecord> = Vec::with_capacity(CHUNK_RECORDS);
        loop {
            chunk.clear();
            while chunk.len() < CHUNK_RECORDS {
                let Some(r) = reader.next() else { break };
                chunk.push(r?);
            }
            if chunk.is_empty() {
                break;
            }
            let out: Vec<(OwnedRecord, u64, u64)> = chunk
                .par_drain(..)
                .map(|mut rec| {
                    let bases_in = rec.seq.len() as u64;
                    process(&mut rec, None, self.cfg)?;
                    let bases_out = rec.seq.len() as u64;
                    Ok((rec, bases_in, bases_out))
                })
                .collect::<Result<Vec<_>>>()?;
            for (rec, bin, bout) in out {
                report.reads_in += 1;
                report.reads_out += 1;
                report.bases_in += bin;
                report.bases_out += bout;
                writer.write_record(&rec.id, &rec.seq, &rec.qual)?;
            }
        }
        writer.finalize()?;
        Ok(report)
    }

    pub fn run_pe(&self, in1: &Path, in2: &Path, out1: &Path, out2: &Path) -> Result<UmiReport> {
        let mut r1 = open_fastq(in1)?;
        let mut r2 = open_fastq(in2)?;
        let mut w1 = ChunkedWriter::create(out1, self.compression)?;
        let mut w2 = ChunkedWriter::create(out2, self.compression)?;
        let mut report = UmiReport {
            mode: Some("PE"),
            input_r1: Some(in1.display().to_string()),
            input_r2: Some(in2.display().to_string()),
            output_r1: Some(out1.display().to_string()),
            output_r2: Some(out2.display().to_string()),
            ..UmiReport::default()
        };
        let mut chunk: Vec<OwnedPair> = Vec::with_capacity(CHUNK_RECORDS);
        let mut done = false;
        while !done {
            chunk.clear();
            while chunk.len() < CHUNK_RECORDS {
                match (r1.next(), r2.next()) {
                    (Some(a), Some(b)) => chunk.push(OwnedPair { r1: a?, r2: b? }),
                    (None, None) => {
                        done = true;
                        break;
                    }
                    _ => {
                        return Err(RsomicsError::InvalidInput(
                            "PE input record counts diverge".into(),
                        ));
                    }
                }
            }
            if chunk.is_empty() {
                break;
            }
            let out: Vec<(OwnedPair, u64, u64)> = chunk
                .par_drain(..)
                .map(|mut p| {
                    let bases_in = (p.r1.seq.len() + p.r2.seq.len()) as u64;
                    process(&mut p.r1, Some(&mut p.r2), self.cfg)?;
                    let bases_out = (p.r1.seq.len() + p.r2.seq.len()) as u64;
                    Ok((p, bases_in, bases_out))
                })
                .collect::<Result<Vec<_>>>()?;
            for (p, bin, bout) in out {
                report.reads_in += 2;
                report.reads_out += 2;
                report.bases_in += bin;
                report.bases_out += bout;
                w1.write_record(&p.r1.id, &p.r1.seq, &p.r1.qual)?;
                w2.write_record(&p.r2.id, &p.r2.seq, &p.r2.qual)?;
            }
        }
        w1.finalize()?;
        w2.finalize()?;
        Ok(report)
    }
}

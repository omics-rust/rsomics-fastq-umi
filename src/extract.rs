use rsomics_common::{Result, RsomicsError};
use rsomics_seqio::OwnedRecord;

use crate::index::{first_index, last_index};
use crate::umi_loc::{UmiConfig, UmiLoc};

// Stamp the UMI tag into the read name before the first space (or append).
pub(crate) fn stamp(id: &[u8], tag: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(id.len() + tag.len());
    if let Some(sp) = id.iter().position(|&b| b == b' ') {
        out.extend_from_slice(&id[..sp]);
        out.extend_from_slice(tag);
        out.extend_from_slice(&id[sp..]);
    } else {
        out.extend_from_slice(id);
        out.extend_from_slice(tag);
    }
    out
}

// fastp trimFront clamps to length()-1 (keeps ≥1 base); trim = min(umi_len+skip, read_len-1).
// Zero-length read has no defined fastp output — fail loud.
pub(crate) fn take_seq_umi(src: &mut OwnedRecord, cfg: &UmiConfig) -> Result<Vec<u8>> {
    let read_len = src.seq.len();
    let umi_len = cfg.len.min(read_len);
    if umi_len == 0 {
        return Err(RsomicsError::InvalidInput(
            "UMI source read has zero length; cannot extract a UMI".into(),
        ));
    }
    let umi = src.seq[..umi_len].to_vec();
    let trim = (umi_len + cfg.skip).min(read_len - 1);
    src.seq.drain(..trim);
    src.qual.drain(..trim);
    Ok(umi)
}

// fastp UmiProcessor::process: if(!umi.empty()) guard means a missing index
// field is a pass-through, not an error.
pub(crate) fn process(
    rec: &mut OwnedRecord,
    mut mate: Option<&mut OwnedRecord>,
    cfg: &UmiConfig,
) -> Result<()> {
    let umi: Vec<u8> = match cfg.loc {
        UmiLoc::Read1 => take_seq_umi(rec, cfg)?,
        UmiLoc::Read2 => {
            let m = mate.as_deref_mut().ok_or_else(|| {
                RsomicsError::ConfigError("--umi_loc read2 requires PE input".into())
            })?;
            take_seq_umi(m, cfg)?
        }
        UmiLoc::Index1 => first_index(&rec.id),
        UmiLoc::Index2 => {
            let m = mate.as_deref().ok_or_else(|| {
                RsomicsError::ConfigError("--umi_loc index2 requires PE input".into())
            })?;
            last_index(&m.id)
        }
        UmiLoc::PerIndex => {
            let mut u = first_index(&rec.id);
            if let Some(m) = mate.as_deref() {
                u.push(b'_');
                u.extend_from_slice(&last_index(&m.id));
            }
            u
        }
        UmiLoc::PerRead => {
            let mut u = take_seq_umi(rec, cfg)?;
            if let Some(m) = mate.as_deref_mut() {
                u.push(b'_');
                u.extend_from_slice(&take_seq_umi(m, cfg)?);
            }
            u
        }
    };
    if umi.is_empty() {
        return Ok(());
    }
    let tag = cfg.tag(&umi);
    rec.id = stamp(&rec.id, &tag);
    if let Some(m) = mate {
        m.id = stamp(&m.id, &tag);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cfg(loc: UmiLoc, len: usize, skip: usize, prefix: &str) -> UmiConfig {
        UmiConfig {
            loc,
            len,
            skip,
            prefix: prefix.as_bytes().to_vec(),
            delimiter: b':',
        }
    }

    fn rec(id: &str, seq: &str, qual: &str) -> OwnedRecord {
        OwnedRecord {
            id: id.as_bytes().to_vec(),
            seq: seq.as_bytes().to_vec(),
            qual: qual.as_bytes().to_vec(),
        }
    }

    #[test]
    fn extract_trims_and_stamps_no_space() {
        let mut r = rec("read1", "ACGTACGTAA", "IIIIIIIIII");
        process(&mut r, None, &cfg(UmiLoc::Read1, 4, 0, "")).unwrap();
        assert_eq!(r.id, b"read1:ACGT");
        assert_eq!(r.seq, b"ACGTAA");
        assert_eq!(r.qual, b"IIIIII");
    }

    #[test]
    fn stamp_inserts_before_first_space() {
        let mut r = rec("read1 1:N:0", "TTTTGGGG", "IIIIIIII");
        process(&mut r, None, &cfg(UmiLoc::Read1, 4, 0, "")).unwrap();
        assert_eq!(r.id, b"read1:TTTT 1:N:0");
        assert_eq!(r.seq, b"GGGG");
    }

    #[test]
    fn skip_removes_extra_bases_after_umi() {
        let mut r = rec("r", "AACCGGTT", "IIIIIIII");
        process(&mut r, None, &cfg(UmiLoc::Read1, 2, 2, "")).unwrap();
        assert_eq!(r.id, b"r:AA");
        assert_eq!(r.seq, b"GGTT");
        assert_eq!(r.qual, b"IIII");
    }

    #[test]
    fn prefix_joined_with_underscore() {
        let mut r = rec("r", "ACGTACGT", "IIIIIIII");
        process(&mut r, None, &cfg(UmiLoc::Read1, 4, 0, "UMI")).unwrap();
        assert_eq!(r.id, b"r:UMI_ACGT");
    }

    // fastp Read::trimFront clamps to length()-1 — a short read always keeps its last base.
    #[test]
    fn umi_len_clamped_keeps_last_base() {
        let mut r = rec("r", "ACG", "III");
        process(&mut r, None, &cfg(UmiLoc::Read1, 8, 0, "")).unwrap();
        assert_eq!(r.id, b"r:ACG");
        assert_eq!(r.seq, b"G");
        assert_eq!(r.qual, b"I");
    }

    #[test]
    fn skip_overrun_keeps_last_base() {
        let mut r = rec("r", "AACCGGTT", "IIIIIIIJ");
        process(&mut r, None, &cfg(UmiLoc::Read1, 4, 20, "")).unwrap();
        assert_eq!(r.id, b"r:AACC");
        assert_eq!(r.seq, b"T");
        assert_eq!(r.qual, b"J");
    }

    #[test]
    fn exact_consume_keeps_last_base() {
        let mut r = rec("r", "AACCGGTT", "IIIIIIIJ");
        process(&mut r, None, &cfg(UmiLoc::Read1, 4, 4, "")).unwrap();
        assert_eq!(r.seq, b"T");
        assert_eq!(r.qual, b"J");
    }

    #[test]
    fn one_base_read_kept_and_stamped() {
        let mut r = rec("r", "A", "I");
        process(&mut r, None, &cfg(UmiLoc::Read1, 1, 0, "")).unwrap();
        assert_eq!(r.id, b"r:A");
        assert_eq!(r.seq, b"A");
        assert_eq!(r.qual, b"I");
    }

    #[test]
    fn empty_source_read_errors() {
        let mut r = rec("r", "", "");
        assert!(process(&mut r, None, &cfg(UmiLoc::Read1, 4, 0, "")).is_err());
    }

    #[test]
    fn index1_stamps_from_header_no_trim() {
        let mut r = rec("read 1:N:0:ACGTAA", "TTTTGGGG", "IIIIIIII");
        process(&mut r, None, &cfg(UmiLoc::Index1, 0, 0, "")).unwrap();
        assert_eq!(r.id, b"read:ACGTAA 1:N:0:ACGTAA");
        assert_eq!(r.seq, b"TTTTGGGG"); // index mode never trims
        assert_eq!(r.qual, b"IIIIIIII");
    }

    #[test]
    fn index_missing_field_is_passthrough_not_error() {
        let mut r = rec("plainname", "ACGT", "IIII");
        process(&mut r, None, &cfg(UmiLoc::Index1, 0, 0, "")).unwrap();
        assert_eq!(r.id, b"plainname"); // empty UMI ⇒ fastp !empty guard ⇒ no stamp
        assert_eq!(r.seq, b"ACGT");
    }

    #[test]
    fn per_index_pe_merges_first_and_last() {
        let mut r1 = rec("p 1:N:0:AAA+CCC", "GGGG", "IIII");
        let mut r2 = rec("p 2:N:0:AAA+CCC", "TTTT", "FFFF");
        process(&mut r1, Some(&mut r2), &cfg(UmiLoc::PerIndex, 0, 0, "")).unwrap();
        assert_eq!(r1.id, b"p:AAA_CCC 1:N:0:AAA+CCC");
        assert_eq!(r2.id, b"p:AAA_CCC 2:N:0:AAA+CCC");
        assert_eq!(r1.seq, b"GGGG");
        assert_eq!(r2.seq, b"TTTT");
    }

    #[test]
    fn per_read_pe_merges_both_seq_umis_and_trims_both() {
        let mut r1 = rec("p 1", "AACCGGGG", "IIIIIIII");
        let mut r2 = rec("p 2", "TTGGCCCC", "FFFFFFFF");
        process(&mut r1, Some(&mut r2), &cfg(UmiLoc::PerRead, 4, 0, "")).unwrap();
        assert_eq!(r1.id, b"p:AACC_TTGG 1");
        assert_eq!(r2.id, b"p:AACC_TTGG 2");
        assert_eq!(r1.seq, b"GGGG");
        assert_eq!(r2.seq, b"CCCC");
    }

    #[test]
    fn read2_and_index2_require_pe() {
        let mut r = rec("r 1:N:0:ACGT", "ACGTACGT", "IIIIIIII");
        assert!(process(&mut r, None, &cfg(UmiLoc::Read2, 4, 0, "")).is_err());
        let mut r2 = rec("r 1:N:0:ACGT", "ACGTACGT", "IIIIIIII");
        assert!(process(&mut r2, None, &cfg(UmiLoc::Index2, 0, 0, "")).is_err());
    }
}

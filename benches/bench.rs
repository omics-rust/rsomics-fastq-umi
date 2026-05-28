use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;
use tempfile::NamedTempFile;

fn bench_fastq_umi(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-fastq-umi");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let fq = manifest.join("tests/golden/se_umi.fastq");
    c.bench_function("rsomics-fastq-umi golden", |b| {
        b.iter(|| {
            let out_file = NamedTempFile::new().unwrap();
            let out = Command::new(black_box(bin))
                .args(["-i", fq.to_str().unwrap(), "-o", out_file.path().to_str().unwrap(), "--umi_loc", "read1", "--umi_len", "8"])
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_fastq_umi);
criterion_main!(benches);

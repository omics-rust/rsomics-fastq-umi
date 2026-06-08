// Hand-computed fastp-0.20.1 UMI golden fixtures: gates correctness where fastp is not installed.
use std::path::PathBuf;
use std::process::Command;

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fastq-umi"))
}

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn se_umi_len6_golden() {
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("o.fq");
    let st = Command::new(ours())
        .args([
            "-i",
            fixture("se_umi.fastq").to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--umi_len",
            "6",
        ])
        .status()
        .unwrap();
    assert!(st.success());

    let expected = "\
@readA:AACCGG 1:N:0:ACGT
TTACGTACGTACGT
+
IIIIIIIIIIIIII
@readB:TTGGCC
AATGCATGCATGCA
+
FFFFFFFFFFFFFF
@readC:GGGGCC desc here
CCAAAATTTTACGT
+
HHHHHHHHHHHHHH
";
    assert_eq!(std::fs::read_to_string(&out).unwrap(), expected);
}

#[test]
fn se_umi_len5_skip2_golden() {
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("o.fq");
    let st = Command::new(ours())
        .args([
            "-i",
            fixture("se_umi.fastq").to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--umi_len",
            "5",
            "--umi_skip",
            "2",
        ])
        .status()
        .unwrap();
    assert!(st.success());

    let expected = "\
@readA:AACCG 1:N:0:ACGT
TACGTACGTACGT
+
IIIIIIIIIIIII
@readB:TTGGC
ATGCATGCATGCA
+
FFFFFFFFFFFFF
@readC:GGGGC desc here
CAAAATTTTACGT
+
HHHHHHHHHHHHH
";
    assert_eq!(std::fs::read_to_string(&out).unwrap(), expected);
}

// fastp 0.20.1 Read::trimFront clamps to length()-1 — reads shorter than umi_len keep their last base.
#[test]
fn se_umi_short_len8_golden() {
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("o.fq");
    let st = Command::new(ours())
        .args([
            "-i",
            fixture("se_umi_short.fastq").to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--umi_len",
            "8",
        ])
        .status()
        .unwrap();
    assert!(st.success());

    let expected = "\
@shortA:AACCG
G
+
I
@shortB:TTG
G
+
F
@shortC:A
A
+
I
";
    assert_eq!(std::fs::read_to_string(&out).unwrap(), expected);
}

// index1: UMI from firstIndex (backward scan to last ':'/'+''), stamped without trimming seq/qual; no index field → pass-through (fastp if(!umi.empty()) guard).
#[test]
fn se_umi_index1_golden() {
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("o.fq");
    let st = Command::new(ours())
        .args([
            "-i",
            fixture("se_umi.fastq").to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--umi_loc",
            "index1",
        ])
        .status()
        .unwrap();
    assert!(st.success());

    let expected = "\
@readA:ACGT 1:N:0:ACGT
AACCGGTTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIII
@readB
TTGGCCAATGCATGCATGCA
+
FFFFFFFFFFFFFFFFFFFF
@readC desc here
GGGGCCCCAAAATTTTACGT
+
HHHHHHHHHHHHHHHHHHHH
";
    assert_eq!(std::fs::read_to_string(&out).unwrap(), expected);
}

// fastp 0.20.1 trimFront throws on zero-length reads; fail loud rather than fabricate a record.
#[test]
fn empty_source_read_errors_cli() {
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("o.fq");
    let st = Command::new(ours())
        .args([
            "-i",
            fixture("se_umi_empty.fastq").to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--umi_len",
            "4",
        ])
        .status()
        .unwrap();
    assert!(
        !st.success(),
        "empty UMI source read must fail loud, not emit a fabricated record"
    );
}

// fastp 0.20.1 rejects PE-only locations (read2/index2/per_index/per_read) without paired input at validation time.
#[test]
fn pe_only_loc_in_se_rejected_cli() {
    for loc in ["per_index", "per_read", "read2", "index2"] {
        let tmp = tempfile::tempdir().unwrap();
        let out = tmp.path().join("o.fq");
        let st = Command::new(ours())
            .args([
                "-i",
                fixture("se_umi.fastq").to_str().unwrap(),
                "-o",
                out.to_str().unwrap(),
                "--umi_loc",
                loc,
                "--umi_len",
                "4",
            ])
            .status()
            .unwrap();
        assert!(
            !st.success(),
            "--umi_loc {loc} in SE must fail loud (fastp rejects it), got success"
        );
    }
}

// Paired-end goldens captured from fastp 0.20.1 (pe_umi_*.r1/.r2). The compat.rs
// byte-equality tests cover the same flags but skip when fastp is off PATH; these
// gate the PE locations in CI where fastp is absent.
fn pe_golden(in1: &str, in2: &str, gold1: &str, gold2: &str, extra: &[&str]) {
    let tmp = tempfile::tempdir().unwrap();
    let o1 = tmp.path().join("o1.fq");
    let o2 = tmp.path().join("o2.fq");
    let mut cmd = Command::new(ours());
    cmd.arg("-i")
        .arg(fixture(in1))
        .arg("-I")
        .arg(fixture(in2))
        .arg("-o")
        .arg(&o1)
        .arg("-O")
        .arg(&o2)
        .args(extra);
    assert!(cmd.status().unwrap().success());
    assert_eq!(
        std::fs::read(&o1).unwrap(),
        std::fs::read(fixture(gold1)).unwrap(),
        "R1 diverges from fastp 0.20.1 golden {gold1}",
    );
    assert_eq!(
        std::fs::read(&o2).unwrap(),
        std::fs::read(fixture(gold2)).unwrap(),
        "R2 diverges from fastp 0.20.1 golden {gold2}",
    );
}

#[test]
fn pe_umi_read1_len8_golden() {
    pe_golden(
        "pe_umi.fastq.r1",
        "pe_umi.fastq.r2",
        "pe_umi_read1_len8.r1",
        "pe_umi_read1_len8.r2",
        &["--umi_len", "8"],
    );
}

#[test]
fn pe_umi_index2_golden() {
    pe_golden(
        "pe_index.fastq.r1",
        "pe_index.fastq.r2",
        "pe_umi_index2.r1",
        "pe_umi_index2.r2",
        &["--umi_loc", "index2"],
    );
}

// per_index merges firstIndex_lastIndex of read1's comment.
#[test]
fn pe_umi_per_index_golden() {
    pe_golden(
        "pe_index.fastq.r1",
        "pe_index.fastq.r2",
        "pe_umi_per_index.r1",
        "pe_umi_per_index.r2",
        &["--umi_loc", "per_index"],
    );
}

// per_read trims umi_len off both mates and stamps the merged umi1_umi2.
#[test]
fn pe_umi_per_read_len6_golden() {
    pe_golden(
        "pe_umi.fastq.r1",
        "pe_umi.fastq.r2",
        "pe_umi_per_read_len6.r1",
        "pe_umi_per_read_len6.r2",
        &["--umi_loc", "per_read", "--umi_len", "6"],
    );
}

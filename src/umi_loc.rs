// fastp 0.20.1 --umi_loc set: Read1/Read2 trim 5' seq; Index1/Index2 use
// read-name index field (no trim); PerIndex = firstIndex_lastIndex;
// PerRead = umi1_umi2, both trimmed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UmiLoc {
    Read1,
    Read2,
    Index1,
    Index2,
    PerIndex,
    PerRead,
}

#[derive(Debug, Clone)]
pub struct UmiConfig {
    pub loc: UmiLoc,
    pub len: usize,      // fastp --umi_len
    pub skip: usize,     // fastp --umi_skip: extra 5' bases removed after the UMI
    pub prefix: Vec<u8>, // fastp --umi_prefix: joined to the UMI with _ when non-empty
    pub delimiter: u8,   // read-name / UMI separator; fastp default :
}

impl UmiConfig {
    pub(crate) fn tag(&self, umi: &[u8]) -> Vec<u8> {
        let mut t = Vec::with_capacity(1 + self.prefix.len() + 1 + umi.len());
        t.push(self.delimiter);
        if !self.prefix.is_empty() {
            t.extend_from_slice(&self.prefix);
            t.push(b'_');
        }
        t.extend_from_slice(umi);
        t
    }
}

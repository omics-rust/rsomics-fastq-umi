// fastp src/read.cpp (MIT): backward scan from len-3, returns the field after
// the last ':' or '+' delimiter in the read name.

// fastp Read::lastIndex: backward from len-3, returns everything after last ':'/'+'
pub(crate) fn last_index(name: &[u8]) -> Vec<u8> {
    let len = name.len();
    if len < 5 {
        return Vec::new();
    }
    for i in (0..=len - 3).rev() {
        if name[i] == b':' || name[i] == b'+' {
            return name[i + 1..len].to_vec();
        }
    }
    Vec::new()
}

// fastp Read::firstIndex: backward from len-3; '+' sets field end to index-1
// (dual-index split), last ':' ends scan; returns the ':'..'+' (or ':'..end) field.
pub(crate) fn first_index(name: &[u8]) -> Vec<u8> {
    let len = name.len();
    if len < 5 {
        return Vec::new();
    }
    let mut end = len;
    for i in (0..=len - 3).rev() {
        if name[i] == b'+' {
            end = i.saturating_sub(1);
        }
        if name[i] == b':' {
            // fastp substr(i+1, end-i) == name[i+1 ..= end] → half-open end+1.
            let stop = (end + 1).min(len);
            let start = (i + 1).min(stop);
            return name[start..stop].to_vec();
        }
    }
    Vec::new()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn index_parsing_matches_fastp_readcpp() {
        // single index: first == last == the trailing field
        assert_eq!(first_index(b"R1 1:N:0:ATCACG"), b"ATCACG");
        assert_eq!(last_index(b"R1 1:N:0:ATCACG"), b"ATCACG");
        // dual index: first = before '+', last = after '+'
        assert_eq!(first_index(b"R1 1:N:0:ATCACG+TGGTCA"), b"ATCACG");
        assert_eq!(last_index(b"R1 1:N:0:ATCACG+TGGTCA"), b"TGGTCA");
        // no delimiter / too short → empty
        assert_eq!(first_index(b"abcd"), b"");
        assert_eq!(last_index(b"readname"), b"");
    }
}

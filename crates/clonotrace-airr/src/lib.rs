//! Strict streaming input and compact records for CT-0.1-native-1.
use flate2::read::MultiGzDecoder;
use serde::Serialize;
use sha2::{Digest, Sha256};
use std::{
    cell::RefCell,
    collections::{BTreeMap, HashSet},
    fs::File,
    io::{self, BufRead, BufReader, Read},
    path::Path,
    rc::Rc,
};

pub type Result<T> = std::result::Result<T, String>;
pub const FIELDS: [&str; 8] = [
    "sequence_id",
    "v_call",
    "j_call",
    "junction",
    "donor_id",
    "sample_id",
    "locus",
    "junction_length",
];
pub const OUTPUT_FIELDS: [&str; 3] = ["clone_id", "clonotrace_status", "clonotrace_reason"];
#[derive(Clone, Debug)]
pub struct Record {
    pub row: usize,
    pub id: String,
    pub scope: (String, String, String),
    pub v: Vec<String>,
    pub j: Vec<String>,
    pub junction: Vec<u8>,
}
#[derive(Clone, Debug, Serialize)]
pub struct Row {
    pub id: String,
    pub reason: Option<&'static str>,
}
#[derive(Clone, Debug)]
pub struct Input {
    pub header: String,
    pub records: Vec<Record>,
    pub rows: Vec<Row>,
    pub sha256: String,
}
#[derive(Clone, Debug)]
pub struct Policy {
    pub changeo: bool,
    pub mappings: BTreeMap<String, String>,
    pub max_rows: usize,
    pub max_junction_length: usize,
}
impl Default for Policy {
    fn default() -> Self {
        Self {
            changeo: false,
            mappings: BTreeMap::new(),
            max_rows: 10_000_000,
            max_junction_length: 10_000,
        }
    }
}
struct HashReader {
    file: File,
    hash: Rc<RefCell<Sha256>>,
}
impl Read for HashReader {
    fn read(&mut self, b: &mut [u8]) -> io::Result<usize> {
        let n = self.file.read(b)?;
        self.hash.borrow_mut().update(&b[..n]);
        Ok(n)
    }
}
pub struct Stream {
    pub reader: Box<dyn BufRead>,
    hash: Rc<RefCell<Sha256>>,
}
impl Stream {
    pub fn digest(&self) -> String {
        format!("{:x}", self.hash.borrow().clone().finalize())
    }
}
pub fn open(path: &Path) -> Result<Stream> {
    let hash = Rc::new(RefCell::new(Sha256::new()));
    let r = HashReader {
        file: File::open(path).map_err(|e| e.to_string())?,
        hash: hash.clone(),
    };
    let reader: Box<dyn BufRead> = if path.extension().is_some_and(|x| x == "gz") {
        Box::new(BufReader::with_capacity(65536, MultiGzDecoder::new(r)))
    } else {
        Box::new(BufReader::with_capacity(65536, r))
    };
    Ok(Stream { reader, hash })
}
pub fn line(reader: &mut dyn BufRead, buffer: &mut Vec<u8>) -> Result<bool> {
    buffer.clear();
    if reader
        .read_until(b'\n', buffer)
        .map_err(|e| e.to_string())?
        == 0
    {
        return Ok(false);
    }
    if buffer.last() == Some(&b'\n') {
        buffer.pop();
        if buffer.last() == Some(&b'\r') {
            buffer.pop();
        }
    }
    if buffer.contains(&b'\r') || buffer.contains(&0) {
        return Err("unsupported embedded CR or NUL".into());
    }
    std::str::from_utf8(buffer).map_err(|_| "input is not UTF-8".to_string())?;
    Ok(true)
}
pub fn calls(value: &str, prefix: &str) -> Option<Vec<String>> {
    let mut genes = Vec::new();
    for call in value.split(',') {
        let (gene, allele) = call
            .split_once('*')
            .map_or((call, None), |(g, a)| (g, Some(a)));
        let suffix = gene.strip_prefix(prefix)?;
        if suffix.is_empty()
            || !suffix
                .bytes()
                .all(|c| c.is_ascii_alphanumeric() || b"/-".contains(&c))
            || allele.is_some_and(|a| a.is_empty() || !a.bytes().all(|c| c.is_ascii_digit()))
        {
            return None;
        }
        genes.push(gene.to_owned());
    }
    Some(genes)
}
pub fn load(path: &Path, policy: &Policy) -> Result<Input> {
    if policy
        .mappings
        .keys()
        .any(|k| !FIELDS.contains(&k.as_str()))
    {
        return Err("unknown canonical mapping field".into());
    }
    let mut stream = open(path)?;
    let mut buffer = Vec::new();
    if !line(&mut *stream.reader, &mut buffer)? {
        return Err("missing header".into());
    }
    let header = String::from_utf8(buffer.clone()).unwrap();
    let cols: Vec<&str> = header.split('\t').collect();
    if cols.iter().any(|s| s.trim().is_empty())
        || cols.iter().copied().collect::<HashSet<_>>().len() != cols.len()
    {
        return Err("blank or duplicate header".into());
    }
    if cols
        .iter()
        .any(|s| OUTPUT_FIELDS.contains(s) || ["CLONE", "cell_id", "CELL"].contains(s))
    {
        return Err(
            "existing output or active single-cell column; rename source annotations explicitly"
                .into(),
        );
    }
    let mut index = BTreeMap::new();
    let mut mapped = HashSet::new();
    for key in FIELDS {
        let name = policy.mappings.get(key).cloned().unwrap_or_else(|| {
            if policy.changeo {
                key.to_uppercase()
            } else {
                key.to_string()
            }
        });
        if !mapped.insert(name.clone()) {
            return Err("two mappings target the same column".into());
        }
        if let Some(i) = cols.iter().position(|&x| x == name) {
            index.insert(key, i);
        } else if key != "junction_length" {
            return Err(format!("missing required column {name}"));
        }
    }
    let mut rows = Vec::new();
    let mut records = Vec::new();
    let mut ids = HashSet::new();
    let mut spans: Vec<std::ops::Range<usize>> = Vec::with_capacity(cols.len());
    while line(&mut *stream.reader, &mut buffer)? {
        if rows.len() >= policy.max_rows {
            return Err("max-rows exceeded".into());
        }
        let text = std::str::from_utf8(&buffer).unwrap();
        spans.clear();
        let mut start = 0;
        for stop in memchr::memchr_iter(b'\t', &buffer) {
            spans.push(start..stop);
            start = stop + 1;
        }
        spans.push(start..buffer.len());
        if spans.len() != cols.len() {
            return Err(format!("ragged row {}", rows.len() + 1));
        }
        let get = |key| &text[spans[index[key]].clone()];
        let id = get("sequence_id");
        if id.trim().is_empty() || !ids.insert(id.to_owned()) {
            return Err("blank or duplicate sequence_id".into());
        }
        if get("donor_id").trim().is_empty() || get("sample_id").trim().is_empty() {
            return Err("blank scope boundary".into());
        }
        let junction = get("junction");
        if junction.len() > policy.max_junction_length {
            return Err("max-junction-length exceeded".into());
        }
        let v = calls(get("v_call"), "IGHV");
        let j = calls(get("j_call"), "IGHJ");
        let reason = if get("locus") != "IGH" {
            Some("outside_IGH_scope")
        } else if junction.is_empty() || !junction.bytes().all(|b| b"ACGT".contains(&b)) {
            Some("empty_or_non_ACGT_junction")
        } else if v.is_none() || j.is_none() {
            Some("missing_or_out_of_profile_VJ_call")
        } else if index
            .get("junction_length")
            .is_some_and(|&i| text[spans[i].clone()] != junction.len().to_string())
        {
            Some("inconsistent_declared_junction_length")
        } else {
            None
        };
        if reason.is_none() {
            records.push(Record {
                row: rows.len(),
                id: id.into(),
                scope: (
                    get("donor_id").into(),
                    get("sample_id").into(),
                    get("locus").into(),
                ),
                v: v.unwrap(),
                j: j.unwrap(),
                junction: junction.as_bytes().to_vec(),
            });
        }
        rows.push(Row {
            id: id.into(),
            reason,
        });
    }
    Ok(Input {
        header,
        records,
        rows,
        sha256: stream.digest(),
    })
}

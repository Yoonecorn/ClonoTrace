//! No-replacement atomic directory transactions and byte-preserving row emission.
use clonotrace_airr::{line, open, Input, Record, Result};
use serde_json::Value;
use std::{
    collections::BTreeSet,
    fs::{self, File},
    io::{BufWriter, Write},
    path::Path,
    time::Instant,
};
fn err(e: impl std::fmt::Display) -> String {
    e.to_string()
}
#[cfg(unix)]
fn publish(from: &Path, to: &Path) -> Result<()> {
    use std::{ffi::CString, os::unix::ffi::OsStrExt};
    let from = CString::new(from.as_os_str().as_bytes()).map_err(err)?;
    let to = CString::new(to.as_os_str().as_bytes()).map_err(err)?;
    #[cfg(target_os = "macos")]
    let status = unsafe {
        libc::renameatx_np(
            libc::AT_FDCWD,
            from.as_ptr(),
            libc::AT_FDCWD,
            to.as_ptr(),
            libc::RENAME_EXCL,
        )
    };
    #[cfg(target_os = "linux")]
    let status = unsafe {
        libc::renameat2(
            libc::AT_FDCWD,
            from.as_ptr(),
            libc::AT_FDCWD,
            to.as_ptr(),
            libc::RENAME_NOREPLACE,
        )
    };
    #[cfg(not(any(target_os = "linux", target_os = "macos")))]
    return Err("atomic exclusive directory publication requires macOS or Linux".into());
    #[cfg(any(target_os = "linux", target_os = "macos"))]
    if status == 0 {
        Ok(())
    } else {
        Err(std::io::Error::last_os_error().to_string())
    }
}
#[cfg(not(unix))]
fn publish(_: &Path, _: &Path) -> Result<()> {
    Err("unsupported atomic publication platform".into())
}
fn writer(path: &Path) -> Result<BufWriter<File>> {
    Ok(BufWriter::with_capacity(
        262144,
        File::create(path).map_err(err)?,
    ))
}
fn finish(mut w: BufWriter<File>) -> Result<()> {
    w.flush().map_err(err)?;
    w.get_ref().sync_all().map_err(err)
}
#[allow(clippy::too_many_arguments)]
pub fn write_with_files(
    input_path: &Path,
    out: &Path,
    input: &Input,
    clones: &[Vec<usize>],
    work: &[Vec<usize>],
    mut manifest: Value,
    extra_files: &std::collections::BTreeMap<String, String>,
) -> Result<()> {
    if fs::symlink_metadata(out).is_ok() {
        return Err("output already exists; replacement is prohibited".into());
    }
    let parent = out
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let parent_dir = File::open(parent).map_err(err)?;
    let temp = tempfile::Builder::new()
        .prefix(".clonotrace-")
        .tempdir_in(parent)
        .map_err(err)?;
    let started = Instant::now();
    let mut labels = vec![None; input.rows.len()];
    for (c, g) in clones.iter().enumerate() {
        for &i in g {
            let slot = &mut labels[input.records[i].row];
            if slot.replace(c + 1).is_some() {
                return Err("duplicate output membership".into());
            }
        }
    }
    if labels.iter().flatten().count() != input.records.len() {
        return Err("retained row conservation failed".into());
    }
    let mut source = open(input_path)?;
    let mut buf = Vec::new();
    if !line(&mut *source.reader, &mut buf)? || buf != input.header.as_bytes() {
        return Err("input changed between passes".into());
    }
    let mut assigned = writer(&temp.path().join("assignments.tsv"))?;
    writeln!(
        assigned,
        "{}\tclone_id\tclonotrace_status\tclonotrace_reason",
        input.header
    )
    .map_err(err)?;
    let mut n = 0;
    while line(&mut *source.reader, &mut buf)? {
        if n >= input.rows.len() {
            return Err("input changed between passes".into());
        }
        assigned.write_all(&buf).map_err(err)?;
        if let Some(label) = labels[n] {
            writeln!(assigned, "\t{label}\tretained\t").map_err(err)?;
        } else {
            writeln!(
                assigned,
                "\t\texcluded\t{}",
                input.rows[n].reason.ok_or("unaccounted row")?
            )
            .map_err(err)?;
        }
        n += 1;
    }
    if n != input.rows.len() || source.digest() != input.sha256 {
        return Err("input changed between passes (SHA-256 mismatch)".into());
    }
    finish(assigned)?;
    let mut summary = writer(&temp.path().join("clones.tsv"))?;
    writeln!(summary,"clone_id\trows\tdistinct_junctions\tv_genes\tj_genes\tjunction_lengths\tdonor_id\tsample_id").map_err(err)?;
    for (c, g) in clones.iter().enumerate() {
        let collect = |f: fn(&Record) -> Vec<String>| {
            g.iter()
                .flat_map(|&i| f(&input.records[i]))
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect::<Vec<_>>()
                .join(",")
        };
        let distinct = g
            .iter()
            .map(|&i| &input.records[i].junction)
            .collect::<BTreeSet<_>>()
            .len();
        let r = &input.records[g[0]];
        writeln!(
            summary,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            c + 1,
            g.len(),
            distinct,
            collect(|r| r.v.clone()),
            collect(|r| r.j.clone()),
            collect(|r| vec![r.junction.len().to_string()]),
            r.scope.0,
            r.scope.1
        )
        .map_err(err)?;
    }
    finish(summary)?;
    let mut diag = writer(&temp.path().join("diagnostics.tsv"))?;
    writeln!(diag, "work_unit\trows\tmethod\tevidence_kind").map_err(err)?;
    for (i, g) in work.iter().enumerate() {
        writeln!(
            diag,
            "{}\t{}\t{}\timplementation_work_unit",
            i + 1,
            g.len(),
            manifest["method"].as_str().unwrap()
        )
        .map_err(err)?;
    }
    finish(diag)?;
    manifest["profile_seconds"]["emit"] = Value::from(started.elapsed().as_secs_f64());
    let mut m = writer(&temp.path().join("manifest.json"))?;
    serde_json::to_writer_pretty(&mut m, &manifest).map_err(err)?;
    writeln!(m).map_err(err)?;
    finish(m)?;
    for (name, contents) in extra_files {
        if !matches!(name.as_str(), "derived.tsv" | "cells.tsv") {
            return Err("unsupported extra output filename".into());
        }
        let mut extra = writer(&temp.path().join(name))?;
        extra.write_all(contents.as_bytes()).map_err(err)?;
        finish(extra)?;
    }
    File::open(temp.path())
        .and_then(|d| d.sync_all())
        .map_err(err)?;
    publish(temp.path(), out)?;
    if let Err(e) = parent_dir.sync_all() {
        eprintln!("clonotrace: output committed; parent directory durability sync failed: {e}");
    }
    Ok(())
}

/// Native wrapper preserves the frozen output behavior.
pub fn write(
    input_path: &Path,
    out: &Path,
    input: &Input,
    clones: &[Vec<usize>],
    work: &[Vec<usize>],
    manifest: Value,
) -> Result<()> {
    write_with_files(
        input_path,
        out,
        input,
        clones,
        work,
        manifest,
        &std::collections::BTreeMap::new(),
    )
}
/// Stream arbitrary supported extension files into the same exclusive transaction.
pub fn transaction(out: &Path, emit: impl FnOnce(&Path) -> Result<()>) -> Result<()> {
    if fs::symlink_metadata(out).is_ok() {
        return Err("output already exists; replacement is prohibited".into());
    }
    let parent = out
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let parent_dir = File::open(parent).map_err(err)?;
    let temp = tempfile::Builder::new()
        .prefix(".clonotrace-")
        .tempdir_in(parent)
        .map_err(err)?;
    emit(temp.path())?;
    File::open(temp.path())
        .and_then(|d| d.sync_all())
        .map_err(err)?;
    publish(temp.path(), out)?;
    if let Err(e) = parent_dir.sync_all() {
        eprintln!("clonotrace: output committed; parent directory durability sync failed: {e}");
    }
    Ok(())
}

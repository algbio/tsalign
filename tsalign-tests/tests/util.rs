use std::{
    env,
    process::{Command, Output, Stdio},
    sync::Mutex,
};

use anyhow::{Result, anyhow};

pub fn run_in_repo_root(args: &str) -> Result<()> {
    let output = run_in_repo_root_output(args)?;
    assert!(output.status.success());
    Ok(())
}

/// Run `cargo run -- <args>` in the repository root and return its output without checking the exit status.
pub fn run_in_repo_root_output(args: &str) -> Result<Output> {
    static DIR_LOCK: Mutex<()> = Mutex::new(());

    {
        let lock = DIR_LOCK.lock().unwrap();
        let current_dir = env::current_dir()?;

        if current_dir.ends_with("tsalign-tests") {
            println!("Current dir: {current_dir:?}");
            let parent_dir = current_dir.parent().ok_or(anyhow!("No parent directory"))?;
            println!("Switching to parent dir: {parent_dir:?}");

            // working directory is this crate, a.k.a. "[...]/template-switch-aligner/tsalign-tests"
            // simulate a call from the repo root by traversing to "../"
            env::set_current_dir(parent_dir)?;
        }
        drop(lock);
    }

    let process = Command::new("cargo")
        .args(["run", "--"].into_iter().chain(args.split_whitespace()))
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()?;
    let output = process.wait_with_output()?;

    println!("{}", output.status);
    println!("stdout:\n{}", String::from_utf8_lossy(&output.stdout));
    println!("stderr:\n{}", String::from_utf8_lossy(&output.stderr));

    Ok(output)
}

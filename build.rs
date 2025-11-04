//! Build script for Rotbond.
//!
//! This script handles Windows-specific resource compilation (e.g., embedding icons and version info)
//! when targeting Windows. It uses runtime checks for the target OS and environment to support
//! cross-compilation and different ABIs like MSVC or GNU (MinGW).
//!
//! # Behavior
//! - On Windows targets, verifies the icon file exists and is readable.
//! - Embeds the icon and metadata into the executable.
//! - For GNU ABI (MinGW), uses `windres` instead of `rc.exe` for resource compilation.
//! - Skips compilation on non-Windows targets.
//!
//! # Dependencies
//! - Requires `winresource` crate for resource handling.
//! - For GNU targets, `windres.exe` must be in PATH (e.g., from MSYS2 MinGW).
//!
//! # Examples
//! Run `cargo build --release` to build with embedded resources.
//!
//! # Notes
//! - Re-runs if `build.rs` or `resources/icon.ico` changes.
//! - Fails gracefully if resources can't be compiled, building without them.

use std::path::Path;

/// Main build function.
///
/// Checks the target OS and environment at runtime to conditionally compile Windows resources.
/// Uses `CARGO_CFG_TARGET_OS` for OS detection and `CARGO_CFG_TARGET_ENV` for ABI (e.g., GNU).
fn main() {
    // Check the target OS at runtime using the environment variable
    let is_windows_target = std::env::var("CARGO_CFG_TARGET_OS")
        .map(|os| os == "windows")
        .unwrap_or(false);

    if is_windows_target {
        println!("cargo:warning=Starting Windows resource compilation...");

        // Re-run if any of these change
        println!("cargo:rerun-if-changed=build.rs");
        println!("cargo:rerun-if-changed=resources/icon.ico");

        // ---- 1. Verify icon exists ----
        let icon_path = "resources/icon.ico";
        println!("cargo:warning=Looking for icon at: {icon_path}");

        if !Path::new(icon_path).exists() {
            eprintln!("cargo:warning=Icon not found at {icon_path}");
            eprintln!("cargo:warning=Building without an icon");
            return;
        }

        // Also check if the file is readable
        match std::fs::metadata(icon_path) {
            Ok(metadata) => {
                println!("cargo:warning=Icon file size: {} bytes", metadata.len());
            }
            Err(e) => {
                eprintln!("cargo:warning=Cannot read icon file metadata: {e}");
                eprintln!("cargo:warning=Building without an icon");
                return;
            }
        }

        println!("cargo:warning=Icon found and readable, proceeding with resource compilation");

        // ---- 2. Build the resource ----
        #[cfg(target_os = "windows")]
        {
            let mut res = winresource::WindowsResource::new();

            // Check target environment for ABI-specific handling
            let target_env = std::env::var("CARGO_CFG_TARGET_ENV").unwrap_or_default();
            if target_env == "gnu" {
                // Use windres for GNU ABI (MinGW) instead of default rc.exe
                res.set_windres_path("windres");
                println!("cargo:warning=Using windres for GNU ABI resource compilation");
            }

            println!("cargo:warning=Setting icon...");
            res.set_icon(icon_path);
            println!("cargo:warning=Icon set successfully");

            // ---- 3. Set version info (0.0.2.0) ----
            // Format: (major << 48) | (minor << 32) | (build << 16) | revision
            let version = (0u64 << 48) | (0u64 << 32) | (2u64 << 16) | (0u64 << 0);
            println!("cargo:warning=Setting version to: 0.0.2.0 (0x{:016x})", version);
            res.set_version_info(winresource::VersionInfo::PRODUCTVERSION, version);
            res.set_version_info(winresource::VersionInfo::FILEVERSION, version);

            // Optional metadata
            res.set("CompanyName", "Le Nhan Pham");
            res.set("FileDescription", "Rotbond: conformer generator");
            res.set("ProductName", "Rotbond");
            res.set("LegalCopyright", "Copyright (c) 2025");

            // ---- 4. Compile ----
            println!("cargo:warning=Compiling resources...");
            match res.compile() {
                Ok(()) => {
                    println!("cargo:warning=Windows resources compiled successfully!");
                    println!("cargo:warning=Icon and metadata should now be embedded in the executable");
                }
                Err(e) => {
                    eprintln!("cargo:warning=Resource compilation failed: {e}");
                    eprintln!("cargo:warning=Make sure the Windows SDK is installed and rc.exe is in PATH (or windres for GNU)");
                    eprintln!("cargo:warning=Building without embedded icon");
                }
            }
        }
        #[cfg(not(target_os = "windows"))]
        {
            println!("cargo:warning=Targeting Windows but building on non-Windows host, skipping resource compilation");
        }
    } else {
        println!("cargo:warning=Skipping Windows resource compilation on non-Windows target");
    }
}
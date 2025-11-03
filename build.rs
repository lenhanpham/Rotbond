// build.rs

#[cfg(target_os = "windows")]
use std::path::Path;

#[cfg(target_os = "windows")]
use winresource;

fn main() {
    #[cfg(target_os = "windows")]
    {
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
        let mut res = winresource::WindowsResource::new();
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
                eprintln!("cargo:warning=Make sure the Windows SDK is installed and rc.exe is in PATH");
                eprintln!("cargo:warning=Building without embedded icon");
            }
        }
    }

    #[cfg(not(target_os = "windows"))]
    {
        println!("cargo:warning=Skipping Windows resource compilation on non-Windows platform");
    }
}
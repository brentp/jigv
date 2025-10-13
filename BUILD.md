# Building jigv

This document provides instructions for building jigv from source on various platforms.

## Automated Builds (GitHub Actions)

The project includes a GitHub Actions workflow that automatically builds binaries for multiple platforms:

- **Linux**: x86_64 and i386
- **macOS**: x86_64 (Intel) and ARM64 (Apple Silicon)
- **Windows**: x86_64

Builds are triggered on:
- Pushes to the master branch
- Pull requests
- Version tags (e.g., `v0.1.11`)
- Manual workflow dispatch

Pre-built binaries from releases can be downloaded from the [Releases page](https://github.com/brentp/jigv/releases).

## Building Locally

### Prerequisites

1. **Nim compiler**: Install Nim 1.6.x or later
   - Linux/macOS: `curl https://nim-lang.org/choosenim/init.sh -sSf | sh`
   - Windows: Download from [nim-lang.org](https://nim-lang.org/install.html)

2. **System dependencies**:
   - **Linux**: `libcurl4-openssl-dev`, `libbz2-dev`, `liblzma-dev`, `zlib1g-dev`
     ```bash
     sudo apt-get install libcurl4-openssl-dev libbz2-dev liblzma-dev zlib1g-dev
     ```
   - **macOS**: Install via Homebrew
     ```bash
     brew install curl-openssl xz bzip2 zlib
     ```
   - **Windows**: Dependencies are typically included with the Nim installation

### Build Steps

1. Clone the repository:
   ```bash
   git clone https://github.com/brentp/jigv.git
   cd jigv
   ```

2. Install Nim dependencies:
   ```bash
   nimble install -y
   ```

3. Build the binary:
   ```bash
   nimble build -d:release
   ```

   The compiled binary will be in the project root as `jigv` (or `jigv.exe` on Windows).

### Platform-Specific Builds

#### Linux i386 (32-bit)

To build a 32-bit binary on a 64-bit Linux system:

1. Install 32-bit development tools:
   ```bash
   sudo dpkg --add-architecture i386
   sudo apt-get update
   sudo apt-get install gcc-multilib g++-multilib libc6-dev-i386
   sudo apt-get install zlib1g-dev:i386 libbz2-dev:i386 liblzma-dev:i386 libcurl4-openssl-dev:i386
   ```

2. Build with 32-bit flags:
   ```bash
   nim c -d:release --cpu:i386 --passC:-m32 --passL:-m32 jigv.nim
   ```

#### macOS ARM64 (Apple Silicon)

To cross-compile for Apple Silicon on an Intel Mac:

```bash
nim c -d:release --os:macosx --cpu:arm64 \
  --passC:"-target arm64-apple-macos11" \
  --passL:"-target arm64-apple-macos11" jigv.nim
```

**Note**: 32-bit (i386) macOS builds are not supported, as Apple discontinued 32-bit app support in macOS 10.15 Catalina.

## Troubleshooting

### Missing Dependencies

If you encounter errors about missing libraries during compilation:
- Make sure all system dependencies are installed
- Check that the Nim package manager (nimble) successfully installed all required packages

### Nim Version Issues

If you encounter compilation errors:
- Ensure you're using a compatible Nim version (1.6.x or later recommended)
- Update nimble packages: `nimble refresh && nimble install -y`

### Cross-Compilation Issues

Cross-compilation can be complex. If you encounter issues:
- Use native compilation on the target platform when possible
- Download pre-built binaries from the [Releases page](https://github.com/brentp/jigv/releases)

## Testing

After building, you can verify the binary works:

```bash
./jigv --help
```

This should display the help message with all available options.

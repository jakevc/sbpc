# Release Process

This document outlines the process for creating new releases of SBPC.

## Versioning

SBPC follows [Semantic Versioning](https://semver.org/):

- **MAJOR** version for incompatible API changes
- **MINOR** version for new functionality in a backward compatible manner
- **PATCH** version for backward compatible bug fixes

## Release Steps

1. Update the version in `Cargo.toml`
2. Update the CHANGELOG.md with the changes in the new version
3. Commit these changes with a message like "Bump version to X.Y.Z"
4. Create and push a tag for the new version:

```bash
git tag -a vX.Y.Z -m "Release vX.Y.Z"
git push origin vX.Y.Z
```

5. The GitHub Actions workflow will automatically:
   - Create a new GitHub Release
   - Build binaries for Linux and macOS (x86_64 and ARM64)
   - Attach the binaries to the release

## Binary Naming Convention

Binaries are named according to the following convention:

```
sbpc-vX.Y.Z-{target}.tar.gz
```

Where `{target}` is one of:
- `x86_64-unknown-linux-gnu` (Linux x86_64)
- `x86_64-apple-darwin` (macOS Intel)
- `aarch64-apple-darwin` (macOS Apple Silicon)

## Installing from Binaries

### Linux (x86_64)

```bash
curl -L https://github.com/jakevc/sbpc/releases/download/vX.Y.Z/sbpc-vX.Y.Z-x86_64-unknown-linux-gnu.tar.gz | tar xz
sudo mv sbpc /usr/local/bin/
```

### macOS (Intel)

```bash
curl -L https://github.com/jakevc/sbpc/releases/download/vX.Y.Z/sbpc-vX.Y.Z-x86_64-apple-darwin.tar.gz | tar xz
sudo mv sbpc /usr/local/bin/
```

### macOS (Apple Silicon)

```bash
curl -L https://github.com/jakevc/sbpc/releases/download/vX.Y.Z/sbpc-vX.Y.Z-aarch64-apple-darwin.tar.gz | tar xz
sudo mv sbpc /usr/local/bin/
```

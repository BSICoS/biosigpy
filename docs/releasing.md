# Releasing Biosigpy

Biosigpy versions are independent from Biosiglib versions. Releases are made
from reviewed commits on `main`; publishing to PyPI is intentionally outside
this procedure.

## Release checklist

1. Update the single version in `src/biosigpy/_version.py`.
2. Move user-visible changes from `Unreleased` into a dated changelog entry and
   add `docs/releases/vX.Y.Z.md`.
3. Confirm that `conformance.json` pins the intended exact Biosiglib commit and
   that Biosiglib has not yet released a different commit for this coordinated change.
4. Run the full Python suite and the pinned Biosiglib manifest validator.
5. Build the documentation with `mkdocs build --strict`.
6. Build both artifacts with `python -m build` from a clean checkout and confirm
   that the source archive contains `conformance.json`, release documentation,
   and the complete test harness.
7. Install the wheel and source archive into separate clean virtual
   environments and verify `biosigpy.__version__`, package metadata, imports,
   and a small algorithm call.
8. Merge the reviewed release PR into `main` and wait for all `main` checks.
9. Create an annotated `vX.Y.Z` tag on that exact merge commit and push it.
10. Verify that the release workflow uploads the wheel, source archive, and
    `SHA256SUMS`, that the checksum entries use release-asset basenames so
    `sha256sum --check SHA256SUMS` works in the download directory, and that
    the GitHub Release uses the committed release notes.

The tag-triggered workflow rejects a tag whose name does not match the package
version and performs both clean installation smoke tests before publishing.

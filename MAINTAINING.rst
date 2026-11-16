Documentation maintenance
=========================

The public website is built from docs/ and mkdocs.yml.
Its GitHub Actions workflow is .github/workflows/docs.yml.

For a local check, install requirements-docs.txt in a separate environment
and run ``python -m mkdocs build --strict``.

Do not commit generated site/ content. On GitHub, configure Pages to use
GitHub Actions; the workflow deploys successful builds from main.

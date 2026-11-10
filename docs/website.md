# Website maintenance

The documentation is prepared for GitHub Pages at
[https://l-gallucci.github.io/Efesto/](https://l-gallucci.github.io/Efesto/).
It uses MkDocs and the Read the Docs theme, following Echidna's documentation setup.

## Activate publishing

After reviewing, committing and pushing the changes yourself:

1. Open the repository **Settings → Pages**.
2. Set **Build and deployment → Source** to **GitHub Actions**.
3. Push to `main`, or run **Actions → Documentation website → Run workflow** on `main`.
4. Open the URL reported by the successful deployment job.

Pull requests build and check the site without publishing. Pushes to `main` build,
check and deploy; a failed build prevents deployment. This workflow publishes an
artifact directly and does not create commits or push a `gh-pages` branch.
The previous separate GitHub Wiki synchronization is retired.

## Maintain sources

| Content | Source |
|---|---|
| User guides and research pages | Markdown in `docs/` |
| Sidebar and theme | `mkdocs.yml` |
| CLI reference | Generated from the public `efesto --help` at build time |
| Library counts and tables | Generated from `hmm_registry.tsv` |
| Logo | Existing `assets/EfestoLogo.png` |
| CSS | `docs/stylesheets/extra.css` |

The build reads maintained sources and a small explicit asset list. It does not
publish ignored local reports, validation raw files or local filesystem paths.
Private developer-only command options are not included in the public CLI reference.

## Local preview

Use a separate documentation environment:

```bash
python3.12 -m venv .venv-docs
source .venv-docs/bin/activate
python -m pip install -r requirements-docs.txt
python -m mkdocs serve
```

For a production check:

```bash
python -m mkdocs build --strict
```

Generated HTML goes to the ignored `site/` directory. Commit the sources and
workflow, not the generated HTML. The documentation environment does not install
bioinformatics dependencies or change the Efesto runtime environment.

Official references: [GitHub Pages custom workflows](https://docs.github.com/en/pages/getting-started-with-github-pages/using-custom-workflows-with-github-pages)
and [MkDocs configuration](https://www.mkdocs.org/user-guide/configuration/).

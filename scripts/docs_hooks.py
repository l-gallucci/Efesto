"""Build public documentation from maintained code and registry sources."""
from collections import Counter
import csv
import os
from pathlib import Path
import posixpath
import re
import subprocess
import sys
from urllib.parse import quote, unquote, urlsplit, urlunsplit

from mkdocs.structure.files import File

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / 'docs'
GENERATED = {'architecture.md': ROOT / 'ARCHITECTURE.md'}
LINK = re.compile(r'(!?\[[^\]]*\]\()([^\s)]+)([^)]*\))')
EXPANSION = (
    'lanM_NF053806', 'lanD_PF41493',
    'wtpA_TIGR03730', 'wtpB_NF040839', 'wtpC_NF040840',
    'tupA_NF041772', 'tupB_NF041773', 'tupC_NF041774',
    'mbnA_TIGR04071', 'mbnB_TIGR04159', 'mbnC_TIGR04160',
)


def enabled(row):
    return row['status'] in {'active', 'experimental'}


def table(headers, rows):
    def cell(value):
        return str(value).replace('|', '\\|').replace('\n', ' ')
    return '\n'.join([
        '| ' + ' | '.join(headers) + ' |',
        '| ' + ' | '.join('---' for _ in headers) + ' |',
        *['| ' + ' | '.join(map(cell, row)) + ' |' for row in rows],
    ])


def catalog_fragments():
    with (ROOT / 'hmm_library/hmm_registry.tsv').open(encoding='utf-8', newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    if not rows or len({r['stem'] for r in rows}) != len(rows):
        raise ValueError('Registry must contain unique model stems')
    for row in rows:
        status = row['status']
        if status not in {'active', 'experimental'} and not status.startswith(('deprecated', 'inactive')):
            raise ValueError('Unrecognized registry status: ' + status)
    counts = Counter(r['status'] for r in rows)
    on = sum(enabled(r) for r in rows)
    summary = (
        f"The registry contains **{len(rows)} profiles**: **{counts['active']} active**, "
        f"**{counts['experimental']} experimental** and **{len(rows) - on} excluded** "
        f"(deprecated/inactive). Ordinary runs search **{on} enabled profiles** before "
        "category selection. Experimental profiles are enabled and require the stated interpretation cautions."
    )
    def count_table(field, title):
        groups = sorted({r[field] for r in rows})
        return table([title, 'Registered', 'Enabled', 'Experimental', 'Excluded'], [
            (group, sum(r[field] == group for r in rows),
             sum(r[field] == group and enabled(r) for r in rows),
             sum(r[field] == group and r['status'] == 'experimental' for r in rows),
             sum(r[field] == group and not enabled(r) for r in rows)) for group in groups
        ])
    def model_table(selected):
        return table(['Model (filename stem)', 'Accession', 'NSEQ', 'Cutoff', 'Status'], [
            (f"`{r['stem']}`", r['acc'] or '—', r['nseq'], r['cutoff'], r['status']) for r in selected
        ])
    by_stem = {r['stem']: r for r in rows}
    return {
        'catalog_summary': summary,
        'category_table': count_table('category', 'Category'),
        'source_table': count_table('source', 'Source label'),
        'fes_table': model_table(sorted((r for r in rows if r['category'] == 'IronSulfur-Assembly'), key=lambda r: r['stem'])),
        'expansion_table': model_table([by_stem[stem] for stem in EXPANSION]),
    }


def on_files(files, config):
    for uri, source in GENERATED.items():
        files.append(File.generated(config, uri, content=source.read_text(encoding='utf-8')))
    files.append(File.generated(config, 'img/logo.png', content=(ROOT / 'assets/EfestoLogo.png').read_bytes()))
    env = os.environ.copy()
    # Public help must stay public even when a maintainer builds in a dev shell.
    env.pop('EFESTO_DEV_MODE', None)
    env['PYTHONPATH'] = str(ROOT / 'src')
    env['COLUMNS'] = '100'
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    help_text = subprocess.run(
        [sys.executable, '-m', 'efesto.cli', '--help'],
        cwd=ROOT, env=env, check=True, capture_output=True, text=True,
    ).stdout.strip()
    if not help_text.startswith('usage:'):
        raise ValueError('Public CLI reference could not be generated')
    content = (
        '# Command reference\n\n'
        'Generated from the public `efesto --help` for this source revision.\n'
        'Use your installed command’s help when working with an older release.\n\n'
        '```text\n' + help_text + '\n```\n'
    )
    files.append(File.generated(config, 'command-reference.md', content=content))
    return files


def on_page_markdown(markdown, page, config, files):
    fragments = catalog_fragments()
    for key, value in fragments.items():
        markdown = markdown.replace('{{' + key + '}}', value)
    if re.search(r'\{\{\w+\}\}', markdown):
        raise ValueError('Unresolved documentation placeholder in ' + page.file.src_uri)
    source = GENERATED.get(page.file.src_uri, DOCS / page.file.src_uri)
    source_to_page = {path.resolve(): uri for uri, path in GENERATED.items()}
    source_to_page[ROOT / 'README.md'] = 'index.md'

    def rewrite(match):
        target = urlsplit(match.group(2))
        if target.scheme or target.netloc or not target.path or target.path.startswith('/'):
            return match.group(0)
        candidate = (source.parent / unquote(target.path)).resolve()
        destination = source_to_page.get(candidate)
        if destination is None and candidate.is_relative_to(DOCS):
            destination = candidate.relative_to(DOCS).as_posix()
        if destination is not None:
            relative = posixpath.relpath(destination, posixpath.dirname(page.file.src_uri) or '.')
            url = urlunsplit(('', '', relative, target.query, target.fragment))
        elif candidate.is_relative_to(ROOT) and candidate.exists():
            url = 'https://github.com/l-gallucci/Efesto/blob/main/' + quote(candidate.relative_to(ROOT).as_posix())
            if target.fragment:
                url += '#' + target.fragment
        else:
            return match.group(0)
        return match.group(1) + url + match.group(3)

    parts = re.split(r'(^```[^\n]*\n.*?^```\s*$)', markdown, flags=re.MULTILINE | re.DOTALL)
    return ''.join(part if part.startswith('```') else LINK.sub(rewrite, part) for part in parts)

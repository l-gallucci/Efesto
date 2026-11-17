# License scope

All Efesto-owned code is licensed under the **GNU Affero General Public License,
version 3 only** (`AGPL-3.0-only`), including the application in `src/efesto/`,
standalone Python/R/shell scripts, setup utilities and tests. Efesto's own project
documentation uses the same license. The full license is in [LICENSE](../LICENSE).

Imported HMMs, reference sequences, source metadata, external software and artwork
are outside this grant. They retain their original providers' terms; see
[third-party notices](../THIRD_PARTY_NOTICES.md), the per-model table in
`../hmm_library/model_licenses.tsv`, and the logo credit in the README.
The six rebuilt HMM contributions explicitly listed in
`../hmm_library/LICENSES.md` use CC BY-SA 4.0 on rights held by Efesto contributors;
their input records retain their original terms.

`GPL-3.0.txt` preserves the upstream GPL version 3 text. `MIT.txt` preserves the
historical Efesto MIT notice: rights already granted for earlier versions remain
available under those earlier grants. It is not a separate current license for
Efesto's standalone tools. Retain applicable copyright and attribution notices.
The CC BY-SA legal text is linked from the HMM license scope.

## Dependencies installed by the user

The wheel and source archive contain Efesto's code and model resources, not the
implementations of Pyrodigal, Pyrodigal-GV, HMMER, pheatmap or plotly. Dependency
declarations and installation instructions do not relicense those projects.
There is no requirement to bundle every dependency's source merely to provide
this source-only workflow. If distributing a combined executable, container or
environment containing dependencies, comply with their applicable licenses and
corresponding-source requirements as well.

The optional R plotting script is a standalone Efesto source file under
AGPLv3. Its separately installed dependencies pheatmap (GPL-2) and plotly (MIT)
retain their original licenses; no package implementation is copied into the
script or shipped by Efesto. Pyrodigal/Pyrodigal-GV are separately installed
GPLv3 Python imports. The current distribution is a source workflow, not a
redistribution of those packages. An import/library call can raise different
questions from invoking an independent executable; these notices do not certify
all possible combined distributions. Anyone distributing such a combination must
assess compatibility and meet the applicable source and notice obligations. See
the [GNU guidance on libraries and interpreters](https://www.gnu.org/licenses/gpl-faq.html#InterpreterIncompat).

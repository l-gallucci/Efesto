"""Package the maintained root HMM library without duplicating source files."""
from pathlib import Path
import csv
import shutil
from setuptools import setup
from setuptools.command.build_py import build_py

ROOT = Path(__file__).resolve().parent


class BuildWithLibrary(build_py):
    def run(self):
        super().run()
        if getattr(self, "editable_mode", False):
            # Editable installations read the maintained repository library.
            return
        destination = Path(self.build_lib) / "efesto" / "data"
        library = ROOT / "hmm_library"
        sources = [p for p in library.iterdir()
                   if p.is_file() and p.suffix in {".tsv", ".txt", ".md"}]
        with (library / "hmm_registry.tsv").open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                path = library / row["hmm_file"]
                if not path.resolve().is_relative_to(library.resolve()) or not path.is_file():
                    raise ValueError("Registry profile missing or outside library: " + row["hmm_file"])
                sources.append(path)
        # Remove stale packaged resources from an earlier build in this tree.
        if destination.exists():
            shutil.rmtree(destination)
        for source in sources:
            target = destination / "hmm_library" / source.relative_to(library)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, target)
        for name in ("LICENSE", "THIRD_PARTY_NOTICES.md", "CITATION.cff"):
            shutil.copyfile(ROOT / name, destination / name)
        examples = destination / "examples"
        examples.mkdir()
        shutil.copyfile(ROOT / "operon_rules.json", examples / "operon_rules.json")


setup(cmdclass={"build_py": BuildWithLibrary})

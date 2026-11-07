#!/usr/bin/env python3
"""Optional maintainer check of installed public headers; core builds need no Python."""
import argparse
import os
from pathlib import Path
import random
import re
import shlex
import subprocess


def main():
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument("include_dir", type=Path, help="installed include/curupixa directory")
  parser.add_argument("--work-dir", type=Path, default=Path("build/header-check"))
  parser.add_argument("--cc", default=os.environ.get("CC", "cc"))
  parser.add_argument("--cxx", default=os.environ.get("CXX", "c++"))
  args = parser.parse_args()
  root = args.include_dir.resolve()
  headers = sorted(root.rglob("*.h"))
  if not headers or not (root / "curupixa.h").is_file():
    parser.error("include_dir must contain the installed curupixa headers")
  graph = {}
  for header in headers:
    children = []
    for name in re.findall(r'^\s*#\s*include\s*"([^"]+)"', header.read_text(), re.M):
      candidates = [header.parent / name, root / name, root / "global" / name]
      child = next((p.resolve() for p in candidates if p.is_file()), None)
      if child is None:
        raise RuntimeError(f"missing local include {name} in {header}")
      children.append(child)
    graph[header] = children
  active, done = set(), set()

  def visit(header, path):
    if header in active:
      raise RuntimeError("include cycle: " + " -> ".join(str(p.relative_to(root)) for p in path + [header]))
    if header in done:
      return
    active.add(header)
    for child in graph[header]:
      visit(child, path + [header])
    active.remove(header)
    done.add(header)

  for header in headers:
    visit(header, [])
  args.work_dir.mkdir(parents=True, exist_ok=True)
  orders = [headers, list(reversed(headers)), random.Random(917).sample(headers, len(headers))]
  checks = 0
  for compiler, standard, suffix in [(args.cc, "c17", "c"), (args.cxx, "c++17", "cpp")]:
    source = args.work_dir / f"public_headers.{suffix}"
    for order in [[h, h] for h in headers] + orders:
      source.write_text("".join(f'#include "{h.relative_to(root)}"\n' for h in order) + "int main(void) { return 0; }\n")
      subprocess.run(shlex.split(compiler) + [f"-std={standard}", "-Wall", "-Wextra", "-Wpedantic", "-Werror",
                     "-fsyntax-only", f"-I{root}", str(source)], check=True)
      checks += 1
  print(f"{len(headers)} installed headers: acyclic; {checks} standalone/repeated/combined C17 and C++17 checks passed")


if __name__ == "__main__":
  main()

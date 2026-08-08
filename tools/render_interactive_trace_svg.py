#!/usr/bin/env python3
"""Render a compiler DOT file to an interactive SVG with hover trace highlighting.

Hovering graph.nodes[i] or any node/edge in one of its dependency traces highlights
that whole trace in red. The script follows the compiler DOT convention that real
edges point from graph/result nodes to AST dependency nodes.

Invisible rank-ordering edges are ignored so layout constraints do not become part
of traces.
"""

from __future__ import annotations

import argparse
from collections import defaultdict, deque
from html import unescape as html_unescape
from pathlib import Path
import re
import subprocess
import sys

QUOTED_ID_RE = re.compile(r'"((?:[^"\\]|\\.)*)"')
SVG_GROUP_RE = re.compile(
    r'(<g id="[^"]+" class="(?:node|edge)")>\n<title>(.*?)</title>',
    re.DOTALL,
)
GRAPH_NODE_RE = re.compile(r"^graph_node_(\d+)$")


def unescape_dot_id(text: str) -> str:
    return text.replace(r'\"', '"').replace(r'\\', '\\')


def extract_real_edges(dot_text: str) -> list[tuple[str, str]]:
    edges: list[tuple[str, str]] = []
    brace_depth = 0
    ignoring_invisible_block_at_depth: int | None = None

    for line in dot_text.splitlines():
        parse_this_line = ignoring_invisible_block_at_depth is None

        if "style=invis" in line:
            parse_this_line = False
            if ignoring_invisible_block_at_depth is None:
                ignoring_invisible_block_at_depth = brace_depth

        if parse_this_line and "->" in line:
            ids = [unescape_dot_id(match) for match in QUOTED_ID_RE.findall(line)]
            for left, right in zip(ids, ids[1:]):
                edges.append((left, right))

        brace_depth += line.count("{") - line.count("}")
        if (
            ignoring_invisible_block_at_depth is not None
            and brace_depth <= ignoring_invisible_block_at_depth
        ):
            ignoring_invisible_block_at_depth = None

    return edges


def trace_from(
    outgoing: dict[str, list[str]], start_node: str
) -> tuple[set[str], set[tuple[str, str]]]:
    seen_nodes: set[str] = set()
    seen_edges: set[tuple[str, str]] = set()
    queue: deque[str] = deque([start_node])

    while queue:
        node = queue.popleft()
        if node in seen_nodes:
            continue
        seen_nodes.add(node)
        for child in outgoing.get(node, ()):
            seen_edges.add((node, child))
            queue.append(child)

    return seen_nodes, seen_edges


def compute_trace_memberships(
    edges: list[tuple[str, str]],
) -> tuple[dict[str, set[int]], dict[tuple[str, str], set[int]]]:
    outgoing: dict[str, list[str]] = defaultdict(list)
    known_nodes: set[str] = set()
    for left, right in edges:
        outgoing[left].append(right)
        known_nodes.add(left)
        known_nodes.add(right)

    graph_indices = sorted(
        int(match.group(1))
        for node in known_nodes
        if (match := GRAPH_NODE_RE.match(node)) is not None
    )

    node_traces: dict[str, set[int]] = defaultdict(set)
    edge_traces: dict[tuple[str, str], set[int]] = defaultdict(set)

    for index in graph_indices:
        start_node = f"graph_node_{index}"
        nodes, trace_edges = trace_from(outgoing, start_node)
        for node in nodes:
            node_traces[node].add(index)
        for edge in trace_edges:
            edge_traces[edge].add(index)

    return node_traces, edge_traces


def svg_title_to_edge(title: str) -> tuple[str, str] | None:
    decoded = html_unescape(title)
    if "->" not in decoded:
        return None
    left, right = decoded.split("->", 1)
    return left, right


def inject_trace_metadata_and_script(
    svg_text: str,
    node_traces: dict[str, set[int]],
    edge_traces: dict[tuple[str, str], set[int]],
) -> str:
    def replace_group(match: re.Match[str]) -> str:
        opening = match.group(1)
        title = html_unescape(match.group(2))

        traces: set[int] = set()
        edge = svg_title_to_edge(match.group(2))
        if edge is not None:
            traces = edge_traces.get(edge, set())
        else:
            traces = node_traces.get(title, set())

        if not traces:
            return match.group(0)

        trace_attr = ",".join(str(trace) for trace in sorted(traces))
        return f'{opening} data-traces="{trace_attr}">\n<title>{match.group(2)}</title>'

    svg_text = SVG_GROUP_RE.sub(replace_group, svg_text)

    script = r'''
<style type="text/css"><![CDATA[
.node.trace-highlight ellipse,
.node.trace-highlight polygon {
  stroke: red !important;
  stroke-width: 4px !important;
  fill: #ffd6d6 !important;
}
.node.trace-highlight text {
  fill: red !important;
  font-weight: bold !important;
}
.edge.trace-highlight path {
  stroke: red !important;
  stroke-width: 3px !important;
}
.edge.trace-highlight polygon {
  stroke: red !important;
  fill: red !important;
}
[data-traces] {
  cursor: pointer;
  pointer-events: all;
}
]]></style>
<script type="application/ecmascript"><![CDATA[
(function () {
  // Robust for standalone SVG files: scripts that are direct children of the
  // root SVG do not reliably have document.currentScript.ownerSVGElement.
  const svg = document.documentElement;
  const byTrace = new Map();
  let lockedTraces = null;

  function addClass(element, className) {
    const classes = new Set((element.getAttribute('class') || '').split(/\s+/).filter(Boolean));
    classes.add(className);
    element.setAttribute('class', Array.from(classes).join(' '));
  }

  function removeClass(element, className) {
    const classes = (element.getAttribute('class') || '').split(/\s+/).filter(Boolean);
    element.setAttribute('class', classes.filter((item) => item !== className).join(' '));
  }

  function tracesOf(element) {
    return (element.getAttribute('data-traces') || '')
      .split(',')
      .filter(Boolean);
  }

  for (const element of svg.querySelectorAll('[data-traces]')) {
    for (const trace of tracesOf(element)) {
      if (!byTrace.has(trace)) {
        byTrace.set(trace, []);
      }
      byTrace.get(trace).push(element);
    }
  }

  function clearHighlight() {
    for (const element of svg.querySelectorAll('.trace-highlight')) {
      removeClass(element, 'trace-highlight');
    }
  }

  function highlight(traces) {
    clearHighlight();
    const highlighted = new Set();
    for (const trace of traces) {
      for (const element of byTrace.get(trace) || []) {
        highlighted.add(element);
      }
    }
    for (const element of highlighted) {
      addClass(element, 'trace-highlight');
    }
  }

  for (const element of svg.querySelectorAll('[data-traces]')) {
    element.addEventListener('mouseover', function () {
      if (lockedTraces === null) {
        highlight(tracesOf(element));
      }
    });
    element.addEventListener('mouseout', function () {
      if (lockedTraces === null) {
        clearHighlight();
      }
    });
    element.addEventListener('click', function (event) {
      event.stopPropagation();
      const traces = tracesOf(element);
      if (lockedTraces !== null && lockedTraces.join(',') === traces.join(',')) {
        lockedTraces = null;
        clearHighlight();
      } else {
        lockedTraces = traces;
        highlight(lockedTraces);
      }
    });
  }

  svg.addEventListener('click', function () {
    lockedTraces = null;
    clearHighlight();
  });
})();
]]></script>
'''

    insert_at = svg_text.rfind("</svg>")
    if insert_at == -1:
        raise ValueError("Graphviz output does not look like SVG: missing </svg>")
    return svg_text[:insert_at] + script + svg_text[insert_at:]


def render_dot_to_svg(dot_path: Path, dot_args: list[str]) -> str:
    command = ["dot", "-Tsvg", *dot_args, str(dot_path)]
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    return result.stdout


def default_output_path(input_path: Path) -> Path:
    return input_path.with_name(f"{input_path.stem}.interactive.svg")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Render DOT to SVG with hover highlighting for graph.nodes[i] traces."
    )
    parser.add_argument("dot_file", type=Path)
    parser.add_argument("--out", "-o", type=Path, help="output SVG path")
    parser.add_argument(
        "--dot-arg",
        action="append",
        default=[],
        help="extra argument passed to dot, e.g. --dot-arg=-Granksep=1.2",
    )
    args = parser.parse_args()

    dot_text = args.dot_file.read_text(encoding="utf-8")
    edges = extract_real_edges(dot_text)
    node_traces, edge_traces = compute_trace_memberships(edges)

    if not node_traces:
        print("error: no graph_node_<i> traces found in DOT file", file=sys.stderr)
        return 1

    svg_text = render_dot_to_svg(args.dot_file, args.dot_arg)
    interactive_svg = inject_trace_metadata_and_script(svg_text, node_traces, edge_traces)

    output_path = args.out or default_output_path(args.dot_file)
    output_path.write_text(interactive_svg, encoding="utf-8")

    trace_count = len({trace for traces in node_traces.values() for trace in traces})
    print(f"wrote {output_path}")
    print(f"embedded hover highlighting for {trace_count} graph.nodes traces")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

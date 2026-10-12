"""Annotation transfer — labels and probabilities. ⛔ This package never touches an embedding.

The seam the design rests on: **nuna owns geometry, bacatlas owns labels.** nuna produces, for every
recipient node, a ranked list of donor nodes with pairwise similarities, and never looks at an
annotation; this package consumes that list, attaches claims, assembles evidence and fits the
accuracy model, and never touches an embedding. The old code sprawled partly because one script
reached into both.
"""

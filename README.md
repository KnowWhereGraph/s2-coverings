# s2-coverings
Tool for creating index-free s2 coverings, at any level

## Background

[S2](http://s2geometry.io/) is a spatial grid system with hierarchy, designed to be easily indexed and queried. Knowledge graphs commonly make use of different geospatial indices for linking spatial data to areas.

At the moment, many graph databases don't have native support for making use of the s2 index system. That's where this tool comes into play.

Rather than relying on geosparql functions (which in turn rely on geosparql support and indices), you can instead pre-materialize the relations between cells and query them though the KnowWhereGraph ontology.

This breaks the reliance on the need for the graph database to support s2 indexing and instead make use of the predicate index from the pre-materialized spatial relations.

## Running

### Docker

Docker should be used to generate the s2 coverings, which can be done with

```bash
docker compose up -d
docker exec -it s2 bash
python3 src/s2.py --level <s2_level>
```

Cell integration can be disabled by adding the `--ni` flag to the command,

```bash
python3 src/s2.py --level <s2_level> --ni
```

A complete list of options can be found by running
```bash
python3 src/s2.py --help
options:
  -h, --help            show this help message and exit
  --level LEVEL         Level at which the s2 cells are generated for
  --format [FORMAT]     The format to write the RDF in. Options are xml, n3, turtle, nt, pretty-xml, trix, trig, nquads, json-ld, hext
  --ni [NI]             When used, s2 integration is disabled
  --compressed [COMPRESSED]
                        use the S2 hierarchy to write a compressed collection of relations at various levels
```
Results will be written to the `output/` folder. The results can then be loaded into your graph database and queried. For more information on querying with the KnowWhereGraph ontology, visit the [docs site](https://knowwheregraph.github.io/#/).

### Locally

Due to the steps involved with installing the s2 libray bindings and different approaches needed for each architecture - running outside of Docker isn't supported. If you're inspired, the Dockerfile has all necessary steps to install the requirements to run the tool.

## Development

### Linting

This repository adheres to the Black tool's default settings. It also makes use of isort for import sorting.  Run these before each pull request. 

**Black and isort are competing for formatting the imports differently**. Run isort _after_ running black.

```commandline
black .
isort .
```

### Tests

Unit tests should be run before each pull request. To run them,

`pytest`

### Contributing

Contributions are welcome! Before working on any new features, please file an issue to make sure the work is in scope with this project. New code should be submitted as a pull request to the `develop` branch.
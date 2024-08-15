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
python3 src/s2.py <s2_level>
```
Results will be written to the `output/` folder. The results can then be loaded into your graph database and queried. For more information on querying with the KnowWhereGraph ontology, visit the [docs site](https://knowwheregraph.github.io/#/).

### Locally

Due to the steps involved with installing the s2 libray bindings and different approaches needed for each architecture - running outside of Docker isn't supported. If you're inspired, the Dockerfile has all necessary steps to install the requirements to run the tool.

## Contributing

Contributions are welcome! Before working on any new features, please file an issue to make sure the work is in scope with this project. New code should be submitted as a pull request to the `develop` branch.
# depytrace

This project provides an algorithm for fast
extraction of high-conductance trees (called *core traces*) 
rooted on designated graph nodes. 
<br>
<br>

**License:** Apache Software License
<br>**Author:** Emmanouil (Manios) Krasanakis
<br>**Dependencies:** `networkx` (required) `pcst_fast` (optional)

# :tent: Roadmap
:x: Method ensemble.<br>
:x: Application for github project understanding.


# :zap: Quickstart
Install the library in your environment,
for example per `pip install depytrace`.
Then, create a `networkx` graph, select a root node and run the snippet:
```python
import depytrace as dp

graph, root = ...
trace = dp.Core()(graph, root)
print(dp.conductance(graph, trace))
```

# :dart: Features
* Near-linear running time with respect 
to the number of edges (doubling the edges, 
approximately doubles running time).
* Provable *1+A* factor approximations of maximum conductance, 
where *A* tend to be small.
* Extensibility to future breakthroughs on Steiner problems.



# :wrench: Customization
Core trace extraction is actually an NP-complete problem. For this
reason, solutions provided by this library are approximate and trade-off
tightness for speed.

In particular, solutions are found 
with an algorithm called *trace contraction*, which internally
relies on iteratively solving an NP-complete problems called *ELOD 
maximization*. In turn, the latter can be translated to rooted Steiner tree 
problems. To accommodate future theoretical breakthroughs,
the library allows setting up custom ELOD solvers, where tighter solvers
translate to tighter theoretical bounds of the trace contraction algorithm.

The dafault ELOD maximizer is a heuristic written in native Python and is
chosen thanks to its cross-platform support. For systems that integrate
the *gcc* compiler (e.g. Linux, Windows with mingw installed) you can 
also use a maximizer provided by the library that adjusts and depends on the 
[pcst_fast](https://github.com/fraenkel-lab/pcst_fast) Steiner tree solver.
To use the related maximizer, just install the *pcst_fast* library 
(if *gcc* is properly set up, this should be as simple as 
`pip install pcst_fast`) and call.

```python
import depytrace as dp

graph, root = ...
trace = dp.Core(dp.cleverRPCST)(graph, root)
print(dp.conductance(graph, trace))
```


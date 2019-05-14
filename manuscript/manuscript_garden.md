---
output:
  pdf_document: default
  html_document: default
---
# Introduction

Succession is a complex ecological process, which includes many interacting abiotic and biotic component. Biotic components of the system create a network of interactions. Functional identity of species involved in interaction network as well as topology the network itself governs many important functions of the ecosystem. Network of interactions affects how the system behaves, based on the functional identity of the species involved. However for the tropical forests topology of these networks and how they depend on the higher trophic levels remain unknown.

Long term debate on the control regime of ecosystems (top down and bottom up) is now considered to be a false dichotomy (Schmitz 2010, and there Hunter and Price 1992). Top-down vs bottom-up together can shape ecosystems. Nevertheless,depending on the species composition and character of interactions between modules in the interaction networks systems can be... Depending on the character of these effects they can cause biomass responses of plant community or remain undetected if total biomass is considered. "Non-consumptive effects have a greater likelihood of propagating strong top-down control over plant abundance and community composition" Schmitz book.... succession in tropical forest can be mainly bottom-up controlled due to very fast growth of early successional plants. Predator removal, experiments however, can cause various consumptive and non-consumptive effects on the herbivores and on plants.

But I would prefer to avoid this distinction. I would rather present it as a phenomenon, whether effects cascade down or not for different groups of organisms. Bottom up and top downmay co-occur but at different scales (Elschot et al. 2017)

Work by Schmitz et al. presented mechanistically how predators can shape the plant community composition and provided some basic understanding of food web dynamics. This provided base for understanding how these processes shape complex ecological networks. There are many studies showing top down direct effects of herbivores and indirect effects of predators on plant communities (citations). However, studies resolving plant-herbivore networks to species level are rather scarce (are they?)

So, what is that I predict? There might not be a top down effect on plants and herbivores. But it might be visible between the components. Envisioning the whole system as a collection of food chains (Schmitz in monograph, p. 17). Is there a method to find these collections?

Top down effects caused by predators can be masked by pooling biomass among multiple food web subsystems (Schmitz, Rooney 2006). Therefore we try to test for the magnitude of those effects after partitioning food wes into its compartments and study top down effect within these compartments.

Top down effects within the compartments? Compartment methods (Zhao, Zhang, Tian, and Xu 2018, compartments based on energy channels. Can I get these channels form allometric relationships of herbivore biomass? I dont think that it could be possible to do for invertebrate predators).

Schmitz suggested that complex networks of interactions can be simplified and analysed within strongly interacting modules (Schmitz Resolving ecosysytem complexity). Therefore to study top down effects we used modules identified using DIRTLPAwb+ for weighted bipartite networks (Weighted modularity Beckett 2016).

Responses of the ecosystem to manipulations can be complex. It has been shown that mortality is higher in the smallest sized groups of insects in the grassland (check that Ovadia and Schmitz 2002). This situation might be different in systems where main predators are birds. This suggests a shift in the size structure of herbivores.

Intra-guild predation enchances biodiversity and functioning in complex ecosystems (Wang, Broose and Gravel 2019, Ecology)

Here we present a study trying to show how food web topology is determined by biotic components of early successional communitites in tropical rain forest.
x. Is the a predictive response of the food web to experimental manipulation of the food chain components? What does change with manipulations of herbivores, predators and pathogenic fungi?
- ANCOVA for vulnerability, generality (specialisation doesnt make sese here because we don't have many species of plants. (vulnerability can tell us if there is reduced variabilty in herbivores feeding on plants. Generality can tell us if under some treatments insects tend to narrow their habitat choice)

x. Are there top-down effects/bottom-up effects on plants and herbivore diversity abundance etc. Not trying to mechanistically explain but rather show presence or absence of these. Can I group insects into modules? Then I can characterise them also by diversity within module and size distributions.
- I think that the plot is ok.

x. We test whether these effects can be approximated looking at the changes in specialization of herbivores within plant herbivore food web. relative to the control. IS there any evidence for the habitat shift for most abundant herbivores. How would omnivory affect prediction from schmitz?
- network figures showing some shift. how this shif can be evaluated?

Explore changes in plant herbivore food web structure subjected to manipulations of top trophic level and introduction of dominant herbivore. in our previous work we shown that predator removal did not have an effect of total productivity of experimental communities however it caused

How this can help us understand succession better? Maybe it will be visible for some species that they shift their feeding domain... they feed on less/more plants...

How to account for the diversity effects?

Again, what could be the purpose of increasing generalist weevil abundance?

What with the presence of omnivores?

Plants vs. herbs? Is one of them reacting in a predictable fashion?
Main  goal is to describe what is changing in this part of the community.


Hypotheses:
- lack of predators can cause

# Materials and methods

# Document options

There are two important files to edit to specify the manuscript information.
First, `authors.yaml` should be self-explanatory; it contains the author names,
email address for the corresponding author, and affiliations. The `infos.yaml`
file is for the manuscript title, keywords, etc. Finally, the `ABSTRACT` file
has the abstract. It can contain markdown formatting.

## Tables

Table legends go on the line after the table itself. To generate a reference to
the table, use `{#tbl:id}` -- then, in the text, you can use `{@tbl:id}` to
refer to the table. For example, the table below is @tbl:id. You can remove the
*table* in front by using `!@tbl:id`, or force it to be capitalized with
`\*tbl:id`.

| Using       |  produces |
|:------------|----------:|
| `@tbl:id`   |   @tbl:id |
| `!@tbl:id`  |  !@tbl:id |
| `\*@tbl:id` | \*@tbl:id |

Table: This is a table, and its identifier is `id` -- we can refer to it using
`{@tbl:id}`. Note that even if the table legend is written below the table
itself, it will appear on top in the compiled document. {#tbl:id}

## Equations

Equations can be referenced using the same syntax as tables, using the `eq`
prefix in place of `tbl`. For example:

$$ y = mx + b $$ {#eq:id}

We can refer to @eq:id in the text.

## Adding references

References go in the `references.json` file, at the root of the project.
References are cited with `@key`, where `key` is the unique identifier of the
reference. Both inline, like @hutc59hsr, and in brackets [@hutc57cr] can be
used.

You can also have footnotes.^[this is a footnote -- it is actually rendered as a side-note in the pre-print format.]

## Figures

Figures can be used with the usual markdown syntax. After the path, you can use
`{#fig:id width=50%}` to specify the width and the reference. See @tbl:id for
how to cite. The code below in the markdown source produces @fig:id.

![This is a figure. Figures can have identifiers, and the width can be changed as well. This legend is a bit long, to show what happens in the preprint mode (it continues in the margin below the limit of the figure).]
(figs/llratio.pdf){#fig:id}
<img src="figs/llratio.pdf" alt="some text"  width="600" height="600">

![Alt](figs/llratio.pdf)
# Other elements

## Code blocks

You can use fenced code blocks to render code:

~~~ javascript
// Update affiliations
var print_affiliations = []
for (var af in affiliations) {
  var afobject = {}
  afobject.id = affiliations[af]
  afobject.text = af
  print_affiliations.push(afobject)
}
~~~

Note that code blocks have line numbers of the left, so this does not interfer
with the line numbers of the text (which are on the right).

## Track changes

You can use `make diff` to create a marked-up pdf document. The git revision can
be specified with the `TAG` variable of `make` (by default, the latest commit).
The other option is `AS`, which can be `draft` or `preprint`, to render the
marked-up version as a draft or as a preprint.

## Editorial marks

[Critic Markup][cm] is rendered:

Don't go around saying{-- to people that--} the world owes you a living. The
world owes you nothing. It was here first. {~~One~>Only one~~} thing is
impossible for God: To find {++any++} sense in any copyright law on the planet.
{==Truth is stranger than fiction==}{>>strange but
true<<}, but it is because Fiction is obliged to stick to possibilities; Truth
isn't.

Note that CriticMarkup is *not* rendered into OpenDocument.

[cm]: http://criticmarkup.com/

## Using with knitr, Weave.jl, ...

Just type `make`. If there is a `Rmd` or `Jmd` document with the same base name,
the makefile will render the markdown document for you.

Note that the extensions *must* be `Rmd` or `Jmd`, with an uppercase first
letter. Of course you will need `knitr` (for `R`) or `Weave.jl` (for `julia`).

Because of the way figures are refered to (using the `@fig:id` syntax), it is
better to generate the figure first, and then call it in the text, using
`fig.show='hide'`. The code below will generate @fig:chunk.


```r
plot(sort(rnorm(200)), type='l')
```

You can then use this figure:

![This is the figure created by the chunck `testfig`, so it is in `figure/testfig-1`. You can use different `dev` in the knitr chunk options, so it is possible to generate pdf or png figures.]

With `knitr`, the `kable` function can create tables. If you add the caption
paragraph immediately below, then these tables can be cited. This is how we
produce @tbl:knit.

| Sepal.Length | Sepal.Width | Petal.Length | Petal.Width | Species |
|-------------:|------------:|-------------:|------------:|:--------|
|          5.1 |         3.5 |          1.4 |         0.2 | setosa  |
|          5.0 |         3.6 |          1.4 |         0.2 | setosa  |
|          5.4 |         3.9 |          1.7 |         0.4 | setosa  |

Table: This is a table, and its identifier is `knit` -- we can refer to it using `{@tbl:knit}`. Note that even if the table legend is written below the table itself, it will appear on top in the compiled document. {#tbl:knit}

# Text example

We posit that four simple rules govern the evolution of networks. First, every
network originally consists of just two species sharing a single interaction;
for example, a plant and its herbivore. Second, a speciation event happens at
the top level (*e.g.* the herbivore) with probability $p$, or at the bottom
level with probability $1-p$. Third, the incipient species starts with all
interactions of its ancestor. Fourth, some of these interactions are lost with
probability $\varepsilon(\lambda, k, c)$, which allows interactions---that are
gained through speciation---to be lost either at a fixed rate $\lambda$ or as a
function of the incipient species' degree $k$. The $c$ parameter modulates this
relationship further by influencing whether high degree of an ancestor
increases, or decreases, the probability of the incipient species losing
interactions. We have used the following formulation for $\varepsilon$:

$$\varepsilon(\lambda, k, c) = \left(1+\left(\frac{1}{\lambda}-1\right)\times c^{k-1}\right)^{-1} \,.  $${#eq:epsilon}

In this formulation, $k$ is the number of interactions of the incipient species,
$\lambda$ is the *basal* rate of interaction loss, and $c$ is a parameter
regulating whether species with more interactions tend to gain or lose
interactions over time. Negative values of $c$ imply that *rich get richer*,
*i.e.* species with more interactions tend to conserve them more over
speciation. The special case of $c = 0$ corresponds to no relationship between
the degree of a species and its probability of losing or retaining an
interaction over speciation. The resulting probability of interaction loss, and
its consequences on degree, is shown in figure. The values of $\varepsilon$
belong to $]0;1[$. Note that, because species are duplicated upon a speciation
event, the network still grows over time. If an incipient species should lose
all of its interactions, then it fails to establish.

These four rules translate directly into steps for the model: pick a level at
random, select a species to duplicate, assess the survival of interactions of
the incipient, and add the incipient to the network. These are performed a fixed
number of time -- we impose an upper limit to the richness at each level, and
when this limit is reached, the incipient species replaces one of the resident
species at random. An equilibrium for the measures of network structure (see
next section) is reached within 1000 timesteps. For all situations, we recorded
the network after 5000 iterations.

## Network measures

### Connectance

Connectance, defined as the ratio of realized interactions on the total number
of potential interactions, is one of the most common descriptor of network
structure. In a bipartite network with $T$ species at the top, and $B$ at the
bottom, having a total of $L$ interactions, it is defined as $Co = L/(T\times
B)$. Connectance has a lower bound, as the network cannot have fewer
interactions that the number of species in its more speciose level -- the
minimal connectance is therefore $c_m = \text{max}(T,B)$. This makes the
connectance of networks of different sizes difficult to compare, especially
since bipartite networks tends to have a low connectance. For this reason, we
used a corrected version of connectance, defined as

$$Co^\star=\frac{L-c_m}{T\times B-c_m} \,.$${#eq:cstar}

This takes values between 0 (the network has the minimal number of interactions)
and 1 (all species are connected), but is robust to variations in species
richness.

# References

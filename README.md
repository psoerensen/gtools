# gtools

Public documentation for statistical genetics tools and teaching materials:
**gsuite**, **gsim**, **gact** and **gteach**.

- [Website](https://psoerensen.github.io/gtools/)
- [gsuite: analysis, tutorials and methods](https://psoerensen.github.io/gtools/gsuite/)
- [gsim: simulation](https://psoerensen.github.io/gsim/)
- [gact: genomic association and annotation](https://psoerensen.github.io/gact/)
- [gteach: genetics, genomics and statistical modelling teaching materials](https://psoerensen.github.io/gteach/)

gsuite is the ecosystem and documentation umbrella for seven public standalone
R packages and native computational libraries. gsim and gact retain their own
public repositories, documentation and release processes. gteach retains its
own teaching repository and website. The seven source
packages are distributed under `site/gsuite/repository`. Mixed-model examples
require the development integration; greml/gsolve R packaging is deferred,
with separate or combined packaging undecided. A root `gsuite` R package
release is not currently planned.

`site/` contains a reviewed, static website snapshot. Website updates and
software releases are separate decisions. Source releases retain their own
version, licence and dependency notices.
The presence of documentation does not change any software licence.

The root homepage is the shared entry point. `site/gsuite/` contains the gsuite
subsite; gsim, gact and gteach retain their existing websites. Older gsuite page URLs
within gtools redirect into the new section.

Pushing changes to `site/` or `.github/workflows/pages.yml` on `main`
automatically publishes the website. The **Publish reviewed website** workflow
can also be run manually on `main`. GitHub Pages must use **GitHub Actions**
as its publishing source. Software release assets do not trigger this workflow.

Site content is maintained in the owning development repositories. Please
report corrections rather than editing generated HTML. No development Git
history is included in this publishing repository.

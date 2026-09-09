# gtools

Public documentation for statistical genetics tools: **gsuite**, **gsim** and
**gact**.

- [Website](https://psoerensen.github.io/gtools/)
- [gsuite: analysis, tutorials and methods](https://psoerensen.github.io/gtools/gsuite/)
- [gsim: simulation](https://psoerensen.github.io/gsim/)
- [gact: genomic association and annotation](https://psoerensen.github.io/gact/)

gsim and gact are already public and retain their own repositories,
documentation and release processes. gsuite is in private development; its
tutorials describe capabilities without providing an installable release.

`site/` contains a reviewed, static website snapshot. Website updates and
software releases are separate decisions. Future approved gsuite releases may
be distributed here, with their own version, licence and dependency notices.
The presence of documentation does not change any software licence.

The root homepage is the shared entry point. `site/gsuite/` contains the gsuite
subsite; gsim and gact retain their existing websites. Older gsuite page URLs
within gtools redirect into the new section.

Pushing changes to `site/` or `.github/workflows/pages.yml` on `main`
automatically publishes the website. The **Publish reviewed website** workflow
can also be run manually on `main`. GitHub Pages must use **GitHub Actions**
as its publishing source. Software release assets do not trigger this workflow.

Site content is maintained in the owning development repositories. Please
report corrections rather than editing generated HTML. No development Git
history is included in this publishing repository.

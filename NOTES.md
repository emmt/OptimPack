# Developer notes for OptimPack

## Publishing a new release

To prepare a new release, follow these steps:

1. Bump version number in [README.md](./README.md) and
   [configure.ac](./configure.ac) and commit the changes.

2. Update the [CHANGES.md](./CHANGES.md) file and commit the changes.

3. Rebuild [configure](./configure) script and make the archive:
   ```
   cd build
   make
   make dist-gzip
   mv optimpack-${VERSION}.tar.gz ../releases/.
   ```

4. Push the changes on GitHub:
   ```
   git push
   ```

5. Open [OptimPack](https://github.com/emmt/OptimPack) page at GitHub and click
   on the [releases](https://github.com/emmt/OptimPack/releases) tab and then
   on the [Draft a new release](https://github.com/emmt/OptimPack/releases/new)
   button.

   - Tag the version.
   - Add a brief description and a detailed description of the most important
     changes.
   - Check/uncheck **This is a pre-release**.
   - Attach the release archive.
   - Click on the green **Publish release** button.

6. In the local git repository, pull the tag with:
   ```
   git pull --all
   ```

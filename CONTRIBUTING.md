# Contributing to SODA

Thank you for helping. SODA is maintained by a small team, so short, clear reports help most.

## Report a problem or ask a question

Open an issue and pick **Bug report** or **Question**. Say which browser, operating system and
data type, and the steps that lead to the problem. Check the README's
[Troubleshooting](README.md#troubleshooting) first: the common failures are there.

**Never paste sample IDs, patient information, passwords, or a link that carries a read pass**
(an Azure `?sv=...&sig=...`, an S3 pre-signed link). Use the synthetic demo
(`python examples/make_demo.py --out demo/raw`) to show a problem when you can.

## Send a change

1. Fork the repository and make a branch.
2. Install with `bash install.sh` or `install.ps1` (see the [README](README.md#quick-start)).
3. Make the change. [ADVANCED.md, Changing the code](docs/ADVANCED.md#7-changing-the-code) says
   which file does what.
4. Run the tests:
   - `python tests/smoke_test.py` must end with **SMOKE TEST PASSED**.
   - For a change to the page, also run the browser tests (Playwright; see
     [REPRODUCE.md, section 7](REPRODUCE.md#7-the-browser-tests)).
   - `python web/SpatialViewer/deid.py check .` must say `[deid] clean`.
5. Update the docs that describe what you changed, and add a line to [CHANGELOG.md](CHANGELOG.md).
6. Open a pull request that says what changed and why, and what you tested.

Keep the docs' style: plain words, short sentences, and say what was tested and what was not.

## Things that must not go into the repository

- The de-identification key (`sample_id_key.csv` or any `*_key.csv`).
- Real data, a prepare cache, or an export of real data.
- Your own `config.js` and logos (they are in `.gitignore`).
- Passwords, storage keys or read-pass links.

## Licence

SODA is released under the Washington University licence in [LICENSE](LICENSE) (see the
[README](README.md#licence)). Ask before you send a large change.

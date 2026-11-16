# Consecutive primes at distance 122,290

**Will Howe, Grant Huffman, Jensen Kohlmeyer, and Liam Kruer**

[Read the proof PDF](consecutive-primes-122290-proof.pdf) · [Download the paper](https://github.com/JENW1N/prime-gap-122290/releases/latest/download/consecutive-primes-122290-proof.pdf)

The three-page mathematical paper includes the full arithmetic certificate as an embedded appendix, the exact starting prime, author contact emails, and an AI acknowledgment.

This package records two proven 8,262-digit consecutive primes at distance 122,290. Every one of the 122,289 interior integers has a verified prime divisor at most 32,467.

The endpoints are p = a + 530520 W and p + 122290. The file `gap_122290_cover_b.json` defines W as the product of its 2,180 prime moduli and a as their simultaneous CRT representative in [0,W). The exact decimal endpoints are in `gap_122290_candidate.json`.

## Verification

Use Python 3.11 or later. Installing `gmpy2` accelerates arithmetic. From this directory run

```sh
python verify_gap_122290.py --left gap_122290_left.ecpp --right gap_122290_right.ecpp --jobs 8
```

The complete rerun on November 16, 2026 checked every interior offset and proved both endpoints in 647 and 659 ECPP steps. The observed result is recorded in `verification-2026-11-16.txt`.

The expanded `interior-divisors.tsv` lists one prime divisor for each interior offset. The compact residue certificate is enough to verify the same coverage.

## Construction

The search optimized admissible residue classes, combined them with the Chinese remainder theorem, presieved the progression parameter, and tested the endpoint pair using CGBN-based GPU arithmetic. The successful search checkpoint records 57,344 processed parameters from 677,077 presieved parameters, 103 left base-two passes, and the paired hit at m = 530520. Endpoint primality was subsequently established by the supplied ECPP certificates.

The `search/` directory preserves the final cover, successful-run parameters and checkpoint, and the original optimization and search scripts. These research scripts retain their original Windows/CUDA environment assumptions. Portable verification is available through the command above. Full discovery timing and the historical 7,995-row reference comparison have not been independently repeated during this audit.

The accompanying paper gives the covering argument and the elliptic-curve primality criterion. Its authors are ordered alphabetically by surname. Prime Gap Records accepted the record on November 16, 2026. The paper has not undergone external mathematical review.

## Context

- [Community table of missing prime gaps](https://primegap-list-project.github.io/lists/missing-gaps/)
- [FrontierShift fs-0428](https://aitonomy.org/projects/open-problems/fs-0428)
- [Erdos Problem 853](https://www.erdosproblems.com/853)
- [Mossinghoff and Trudgian, Frangible prime gaps](https://arxiv.org/abs/2610.34525)
- [NVIDIA CGBN](https://github.com/NVlabs/CGBN)

The theorem establishes a specific occurrence. The least occurrence, the prime index, and the asymptotic questions in Erdos Problem 853 remain undetermined.

## Attribution and status

Authors in alphabetical order by surname are Will Howe, Grant Huffman, Jensen Kohlmeyer, and Liam Kruer. The Prime Gap Records team identifier `JKWHLKGH` uses each author's initials to fit its eight-character discoverer field. The work was developed in an AI-assisted Codex workflow. Prime Gap Records lists the result as `122290 6.4287 JKWHLKGH`, with classification `C?D`. Its server reported `Verified! (10.8s)`. See [the live status page](https://primegaps.cloudygo.com/status) and `prime-gap-records-receipt.json`. FrontierShift review is pending.

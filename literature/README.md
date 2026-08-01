# Campaign literature

Reference papers for the DNS validation campaign (and successors). PDFs in
this directory are **untracked** (copyright — see the `.gitignore` rule);
this index is the committed record of what should be here. Full citations
with campaign roles: `dns-validation-campaign-plan.md` §12. EL-era PDFs
remain at repo root until the dashboard's literature-path check is updated.

Naming convention: `surname[_surname...]_year.pdf`, lowercase.

| File | Citation | Campaign role |
|---|---|---|
| `wan_turek_2006.pdf` | Wan & Turek, IJNMF 51 (2006) 531–566 | FBM method lineage (force integral, collision model) |
| `munster_mierka_turek_2012.pdf` | Münster, Mierka & Turek, IJNMF 69 (2012) 294–313 | 3D FEM-FBM lineage of this codebase |
| `hasimoto_1959.pdf` | Hasimoto, JFM 5 (1959) 317–328 | D1 array-drag expansion |
| `zick_homsy_1982.pdf` | Zick & Homsy, JFM 115 (1982) 13–26 | D1 array drag, general φ |
| `uhlmann_dusek_2014.pdf` | Uhlmann & Dušek, IJMF 59 (2014) 221–243 | D1.4 free-fall benchmark + resolution requirements |
| `brenner_1961.pdf` | Brenner, Chem. Eng. Sci. 16 (1961) 242–251 | D2 sphere–wall approach |
| `cooley_oneill_1969.pdf` | Cooley & O'Neill, Mathematika 16 (1969) 37–49 | D2 near-contact asymptotics |
| `jeffrey_onishi_1984.pdf` | Jeffrey & Onishi, JFM 139 (1984) 261–290 | D2 two-sphere resistance functions |
| `fortes_joseph_lundgren_1987.pdf` | Fortes, Joseph & Lundgren, JFM 177 (1987) 467–483 | D2 DKT experiment |
| `glowinski_pan_hesla_joseph_periaux_2001.pdf` | Glowinski et al., JCP 169 (2001) 363–426 | D2 DKT computation |
| `apte_martin_patankar_2009.pdf` | Apte, Martin & Patankar, JCP 228 (2009) 2712–2738 | D2 DKT candidate primary reference |
| `breugem_2012.pdf` | Breugem, JCP 231 (2012) 4469–4498 | D2 DKT alternative reference |
| `beestra_2007.pdf` | Beetstra, van der Hoef & Kuipers, AIChE J. 53 (2007) 489–501 | D3.1 random-array drag correlation (filename keeps provider's spelling) |
| `tenneti_2011.pdf` | Tenneti, Garg & Subramaniam, IJMF 37 (2011) 1072–1092 | D3.1 second drag correlation |
| `richardson_1954.pdf` | Richardson & Zaki, Trans. Inst. Chem. Eng. 32 (1954) 35–53 | D3.2 hindered-settling exponent band |
| `ergun_1952.pdf` | Ergun, Chem. Eng. Prog. 48 (1952) 89–94 | D3.3 fixed-bed Δp / u_mf anchor |
| `jeffrey_1922.pdf` | Jeffery, Proc. R. Soc. Lond. A 102 (1922) 161–179 (title page verified) | D4.2 orbit period (filename keeps provider's spelling; author is G.B. Jeffery) |
| `happel_brenner_1973.pdf` | Happel & Brenner, *Low Reynolds Number Hydrodynamics*, Noordhoff (full book) | D4.3 spheroid drag formulas |
| `ardekani_2016.pdf` | Ardekani, Costa, Breugem & Brandt, IJMF 87 (2016) 16–34 (arXiv 1602.05769 preprint; title page verified) | D4-ORIENT primary; also covers **spheroid DKT** — bridges D2.3/D4 |
| `causin_2005.pdf` | Causin, Gerbeau & Nobile, CMAME 194 (2005) 4506–4527 (title page verified) | Theory for the D1 dt **stability floor** finding (datasheet `dt_stability`): added-mass instability of loosely-coupled partitioned schemes, aggravated by decreasing dt |

Still wanted: **Di Felice 1994** (D5 closure under test) — the only
outstanding paper. Note: `happel_brenner_1073.pdf` is a byte-identical
duplicate of `happel_brenner_1973.pdf` (same md5) and can be deleted.

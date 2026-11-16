from pathlib import Path
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, PageBreak, Table, TableStyle, KeepTogether
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from xml.sax.saxutils import escape

out=Path('outputs')
styles=getSampleStyleSheet()
styles.add(ParagraphStyle(name='PaperTitle',fontName='Times-Bold',fontSize=20,leading=23,spaceAfter=15,alignment=TA_CENTER))
styles.add(ParagraphStyle(name='PaperBody',fontName='Times-Roman',fontSize=11,leading=14,spaceAfter=9))
styles.add(ParagraphStyle(name='PaperSmall',fontName='Times-Roman',fontSize=9,leading=11,spaceAfter=6))
styles.add(ParagraphStyle(name='Equation',fontName='Courier',fontSize=10,leading=14,spaceBefore=6,spaceAfter=12,leftIndent=16))
styles['Heading1'].fontName='Times-Bold';styles['Heading1'].fontSize=14;styles['Heading1'].leading=17
story=[]
def p(t,style='PaperBody'): story.append(Paragraph(t,styles[style]))
def h(t): p(t,'Heading1')
def eq(t):p(escape(t),'Equation')
def page():story.append(PageBreak())
p('A certified consecutive-prime gap of length 122,290','PaperTitle')
p('Jensen Kohlmeyer, Will Howe, Liam Kruer, and Grant Huffman<br/>Research note draft | 16 November 2026','PaperSmall')
h('Abstract')
p('We exhibit consecutive primes p and p + 122,290, each with 8,262 decimal digits. A finite covering system using 2,180 distinct prime moduli, with largest modulus 32,467, supplies a proper prime divisor for every integer strictly between the endpoints. Elliptic-curve primality certificates prove both endpoints prime. The construction combines a search for economical residue covers, a Chinese remainder theorem progression, sieving of its parameter, and GPU-assisted endpoint testing. The supplied verifier reconstructs the progression and checks the entire certificate. This gives an occurrence of the length identified as missing in the Prime Gap List and in FrontierShift task fs-0428. Public record acceptance remains pending.')
h('1. The result')
p('<b>Theorem 1.</b> Let Q be the 2,180 primes listed in <font face="Courier">gap_122290_cover_b.json</font>, and let b_q denote the accompanying residue for q. Define W as their product and a as the unique integer in [0,W) satisfying a = b_q modulo q for every q in Q. Then')
eq('p = a + 530520 W,    q = p + 122290')
p('are consecutive primes. Both p and q have 8,262 decimal digits. Every integer strictly between them has a prime divisor at most 32,467. In particular this is a frangible prime gap in the terminology of Mossinghoff and Trudgian [3].')
p('The JSON endpoint file records both integers in decimal. W has 8,256 decimal digits. The gap merit 122290/log(p) is approximately 6.428675292687395, where log denotes the natural logarithm.')
p('The data, certificates, verification program, and expanded list of interior divisors form the finite specification of the theorem. Hashes identifying the endpoint strings and certificate files appear in Section 5.')
h('2. The record question')
p('The Prime Gap List tracks known occurrences of specified consecutive-prime gap lengths [1]. FrontierShift fs-0428 asks for an occurrence of exactly 122,290, with full interior verification [2]. On 16 November 2026, the live merits file contained 126,372 entries and no entry of size 122,290. This establishes the status of that catalog at the time of our check. Submission and community review are separate steps.')
p('The wider question of how the smallest missing even gap changes with the range of prime indices appears as Erdos Problem 853 [4]. The present theorem supplies one additional gap length. The asymptotic questions recorded there remain open.')
page()
h('3. Constructing the arithmetic progression')
p('Write h = 122,290 and H = h/2 = 61,145. We require odd endpoints. The congruence a = 1 modulo 2 consequently covers every odd offset j, because a+j is then even. For even offsets j=2t, with 1 &lt;= t &lt; H, choose odd primes q and residues r_q such that')
eq('t = r_q (mod q),    r_q != 0,    r_q != H (mod q).')
p('The first condition identifies sites covered by q. The two exclusions preserve the sites t=0 and t=H corresponding to the endpoints. Convert the selected half-interval residues to endpoint residues by setting b_q = -2r_q modulo q. Then a+2t is divisible by q at every covered site. Distinct prime moduli allow the Chinese remainder theorem to combine the assignments.')
p('The saved program <font face="Courier">optimize_gap_mask.py</font> constructs candidate covers using randomized coordinate updates. It maintains the number of primes covering each site. When updating one prime, it temporarily removes its residue class and chooses an admissible replacement that covers as many currently uncovered sites as possible. Randomized order, tie-breaking, and removal of small collections of classes allow interacting assignments to change together.')
p('The saved repair program groups remaining sites by admissible larger-prime residue classes. It greedily favors groups with high coverage relative to log(q), then matches the remaining sites to unused prime moduli. The improvement program uses bounded CP-SAT subproblems to replace subsets of classes while minimizing the sum of log(q). This objective corresponds to log(W), so smaller values reduce the size of every endpoint subsequently tested. No global optimality of the final cover is claimed.')
p('The final cover is saved as <font face="Courier">improved_122290.json</font>. Converting its half-interval assignments gives exactly the residue list in <font face="Courier">gap_122290_cover_b.json</font>. The final certificate, including all its coverage checks, is sufficient to reconstruct the result independently of the randomized search that produced it.')
h('3.1 A finite covering lemma')
p('<b>Lemma 2.</b> Suppose Q is a finite set of distinct primes, W is their product, and a satisfies a = b_q modulo q for q in Q. Suppose every offset 1 &lt;= j &lt; h has some q in Q with b_q+j = 0 modulo q. If p=a+mW exceeds every q in Q and both p and p+h are prime, then they are consecutive primes.')
p('<b>Proof.</b> For each interior j, the chosen q divides a+j and W, hence divides p+j. Also 1 &lt; q &lt; p+j. Thus p+j is composite. No integer strictly between the prime endpoints is prime, establishing consecutiveness. If all q &lt; h, the gap is frangible. This last condition holds for the certificate here, since 32,467 &lt; 122,290. &#9633;')
page()
h('4. Searching for prime endpoints')
p('After the cover is fixed, every parameter m produces the same forced-composite interior in the interval from a+mW to a+mW+h. The remaining task is to find a parameter giving two prime endpoints. For a sieving prime ell outside Q, W is invertible modulo ell. The two excluded parameter classes are')
eq('m = -a/W (mod ell),    m = -(a+h)/W (mod ell).')
p('The preparation code marks these classes before performing large-integer tests. The successful recorded run used 300,001 &lt;= m &lt;= 3,100,000 and sieving primes through 500,000,000. Its metadata report 677,077 eligible parameters. The saved GPU checkpoint records 57,344 processed parameters, 103 left endpoints passing the base-two test, and one double-passing parameter, m=530,520.')
p('The GPU code uses NVIDIA CGBN multiple-precision arithmetic and Montgomery modular squaring [5]. It evaluates the base-two Fermat condition for left endpoints, then for corresponding right endpoints that remain. The successful candidate was independently checked with GMP probable-prime tests and subsequently certified by ECPP. Fermat and probable-prime testing serve as search filters. The final theorem is established by the certificates described in Section 5.')
p('The retained project report identifies an RTX 5090 and records agreement with PFGW on a separate 7,995-row reference shard. The current publication package preserves the successful-run parameters and checkpoint. The reference comparison is reported as provenance from that earlier run, and was not repeated during the 16 November audit. Likewise, elapsed search time is omitted because a complete timing record has not been established for every preparation and search phase.')
p('This procedure contains a computational search over both residue covers and progression parameters. The congruence design targets the exact desired gap and makes interior compositeness automatic for every eligible parameter. Endpoint primality is then located by enumerating and filtering candidates.')
h('4.1 Reproducibility scope')
p('The certificate verifier uses the Python standard library and optionally gmpy2 for accelerated arithmetic. It accepts the saved residue list and endpoint certificates directly. Reproducing discovery additionally requires NumPy, Numba, SymPy, OR-Tools, CUDA bindings, and the CGBN headers and setup expected by the supplied search scripts. The latter scripts preserve the original research environment assumptions; the portable verification command is the primary reproducibility entry point.')
p('The successful-run checkpoint is a record of a completed search segment. It does not establish that 530,520 is the first parameter yielding prime endpoints, and the result does not identify the globally least occurrence of this gap size.')
page()
h('5. Proving and checking the certificate')
p('Interior verification first proves each listed modulus prime by trial division, verifies distinctness and the endpoint exclusions, and reconstructs a by successive Chinese remainder theorem updates. It checks all 122,289 interior offsets against the cover. The expanded TSV file additionally supplies a prime divisor for each offset. Decimal endpoints must match the reconstruction and their recorded SHA-256 digests.')
p('The left ECPP chain has 647 steps and terminates at 77,611,968,613. The right chain has 659 steps and terminates at 78,278,019,539. The terminal primes are below 2^64 and are checked with a deterministic Miller-Rabin base set. The verifier then certifies every chain step using exact integer arithmetic and elliptic-curve group operations.')
p('For a step with modulus n, point P, cofactor s, and recursively certified prime r, it reconstructs the curve, verifies nonsingularity modulo each prime divisor of n through a greatest-common-divisor check, and computes R=[s]P. Every affine addition denominator must be a unit modulo n. It requires R to be finite and [r]R to be the identity, and checks r &gt; (n^(1/4)+1)^2 using integer arithmetic. If n were composite, reduction modulo a prime divisor ell &lt;= sqrt(n) would produce a point of order r. Hasse\'s bound gives r &lt;= (sqrt(ell)+1)^2 &lt;= (n^(1/4)+1)^2, a contradiction. Thus n is prime. This is the standard ECPP certification principle [6].')
p('On 16 November 2026, the complete supplied Python verifier was rerun with eight worker processes and exited successfully after proving both endpoints and checking the full interior. The preceding verification log is included in the package. A separate PARI/GP run started during that audit was stopped after the Python verification completed; no completed result from that additional run is asserted here.')
eq('python verify_gap_122290.py --left gap_122290_left.ecpp\n       --right gap_122290_right.ecpp --jobs 8')
p('SHA-256 digests of ASCII endpoint decimal strings, without a newline','PaperSmall')
for s in ['p  2d7a21344e3dd26dfb667033b31d7f40a1459b266b3c837ded335fd3e0a27d05','q  b450c29e7ac6bbe29083f66901baeccdd30e24fb3c2eb5ac82afc02ad1982a47']:
 p(s[:35]+'<br/>'+s[35:],'PaperSmall')
p('SHA-256 digests of certificate files','PaperSmall')
for name,s in [('Left','402af7bca1189862e49075219d2609aec7f1e726a25b71a0f74ae16e50b70213'),('Right','6c3aa8c159124ae6841271026c61723cd2bde0f2301ce2dd25d43f929e9db91f')]:
 p(name+'<br/>'+s[:32]+'<br/>'+s[32:],'PaperSmall')
page()
h('6. Interpretation and further work')
p('The main result is a new explicitly certified occurrence of one prescribed gap length. The method consists of finite covering optimization, arithmetic-progression sieving, endpoint search, and independent certificate verification. Chinese remainder constructions and residue-cover descriptions of frangible gaps have established precedents [3]. A claim of general methodological novelty would require a broader comparison with earlier exact-gap and covering searches.')
p('A useful next study would compare complete covers with partial covers under fixed computational budgets, measuring both the cost of endpoint tests and the cost of certifying uncovered interior sites. The relevant quantities include log(W), the proportion of parameters surviving the sieve, left and paired endpoint yield, certificate-generation cost, and full verification cost. This note supplies one verified instance and the search artifacts for examining it.')
h('Computational assistance and status')
p('The construction and this draft were developed in an AI-assisted Codex workflow. Search and verification programs, finite witnesses, and audit results are supplied for independent scrutiny. This draft has not undergone external mathematical review. The public catalog status and FrontierShift acceptance should be updated from their eventual recorded responses.')
h('References')
refs=[
('[1]','Prime Gap List Project. Missing prime gaps; data and submission instructions. Accessed 16 November 2026.','https://primegap-list-project.github.io/lists/missing-gaps/'),
('[2]','AItonomy Foundation. FrontierShift fs-0428, A prime gap of exactly 122290. Screened 5 August 2026.','https://aitonomy.org/projects/open-problems/fs-0428'),
('[3]','M. J. Mossinghoff and T. S. Trudgian. Frangible prime gaps. arXiv:2610.34525, 2026.','https://arxiv.org/abs/2610.34525'),
('[4]','T. F. Bloom. Erdos Problem 853. Accessed 16 November 2026.','https://www.erdosproblems.com/853'),
('[5]','NVIDIA. CGBN, Cooperative Groups Big Numbers. Source repository.','https://github.com/NVlabs/CGBN'),
('[6]','FLINT documentation. Elliptic curve primality proving and certificate verification.','https://flintlib.org/doc/ecpp.html'),
]
for n,t,u in refs: p(f'{n} {escape(t)}<br/><link href="{u}" color="#224488">{escape(u)}</link>','PaperSmall')
def footer(canvas,doc):
 canvas.setFont('Times-Roman',9);canvas.drawString(52,28,'Certified gap 122,290 | Research note draft');canvas.drawRightString(560,28,str(doc.page))
doc=SimpleDocTemplate(str(out/'gap-122290-paper-draft.pdf'),pagesize=(612,792),rightMargin=52,leftMargin=52,topMargin=45,bottomMargin=45)
doc.build(story,onFirstPage=footer,onLaterPages=footer)
print(out/'gap-122290-paper-draft.pdf')

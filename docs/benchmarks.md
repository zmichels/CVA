# Bounded CVA timing observations

All measurements below use the same inputs and options within each comparison. They are wall-clock observations on one macOS ARM64 workstation running MATLAB R2025b with MTEX 6.1.1, not general throughput guarantees. The exact mean-referenced method is compared with itself; the scientifically distinct center-reference shortcut is excluded from the speedup claims.

| Case | Legacy/reference times (s) | Maintained batched times (s) | Median ratio |
| --- | --- | --- | ---: |
| Grain CVA, 24 selected grains across three phases, no KDE or stored tensor, one MATLAB single-thread run | 0.159, 0.160, 0.103 | 0.110, 0.130, 0.076 | 1.44× |
| Mean-referenced grid CVA, 324 kernels in a 20-by-20 forsterite crop, including identical 0.25°/10° KDE, three alternating root/new calls | 3.087, 1.996, 1.144 | 1.008, 0.821, 0.748 | 2.43× |

The grid experiment warmed up substantially across repetitions; both functions ran in alternating order and the median is more representative than the first ratio. Its process printed all timing and parity observations, then MATLAB crashed during native batch shutdown after the density work, so the result has a process-exit qualification. The clean no-density grid test subset exits 0. The grain timings are a small sample and include per-grain MTEX mean calls in both methods; most batched savings come from fewer per-pixel log and tensor operations.

The earlier ODA research prototype measured approximately 7.4× on a different 80-grain case by reusing grain-object means. The released maintained grain path deliberately computes the same per-grain MTEX means as the root legacy method; the 7.4× figure is **not** its benchmark.

MTEX for Python calls retain symmetry-aware means and specimen-frame logs. This release establishes numerical parity first and makes no Python-versus-MATLAB speed claim. Benchmark a larger representative crop with identical phase and grain selections before using either implementation's wall time to size a production workload.

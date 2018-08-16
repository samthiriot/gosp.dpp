
#source('functions.R')

#source('data1.R')

data(dwellings_households)

# compute some statistics on the input data and parameters
prepared <- matching.prepare(dwellings_households$sample.A, dwellings_households$sample.B, dwellings_households$pdi, dwellings_households$pdj, dwellings_households$pij)

# define how to distribute biases based on the input data and control parameters

#fonctionnel:
#disc <- matching.solve(case.prepared, nA=50000,nB=40000, nu.A=0, phi.A=0, delta.A=1, nu.B=0, phi.B=0, delta.B=0, gamma=1)
#disc <- matching.solve(case.prepared, nA=50000,nB=40000, nu.A=0, phi.A=0, delta.A=1, gamma=1, delta.B=1, phi.B=0, nu.B=1)
#disc <- matching.solve(case.prepared, nA=50000, nB=40000, nu.A=0, phi.A=0, delta.A=0, gamma=1, delta.B=1, phi.B=1, nu.B=1)
solved <- matching.solve(prepared, nA=50000, nB=40000, nu.A=1, phi.A=1, delta.A=1, gamma=0, delta.B=1, phi.B=1, nu.B=1)


# refusé à raison:
# si on ne laisse pas assez de flexibilité (on retire flex sur delta A)
#disc <- matching.solve(case.prepared, nA=50000,nB=40000, nu.A=0, phi.A=0, delta.A=0, nu.B=0, phi.B=0, delta.B=0, gamma=1)
# si on ne laisse pas assez de flex sur p
#disc <- matching.solve(case.prepared, nA=50000,nB=40000, nu.A=0, phi.A=1, delta.A=0, nu.B=0, phi.B=0, delta.B=0, gamma=0)

# TODO cas dans lesquels il reste un doute; on peut répartir le risque en fonction des poids

# disc <- matching.solve(case.prepared, nA=50000, nB=40000, nu.A=1, phi.A=0, delta.A=1, gamma=0, delta.B=1, phi.B=1, nu.B=1)
print(disc)

# generate the population
generated <- matching.generate(solved, dwellings_households$sample.A, dwellings_households$sample.B)

plot(generated, dwellings_households$sample.A$sample, dwellings_households$sample.B$sample)

print(case)


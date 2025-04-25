commutator(A, B) = A * B - B * A

anticommutator(A, B) = A * B + B * A

vectomat( vec ) = mapreduce( permutedims, vcat, vec)
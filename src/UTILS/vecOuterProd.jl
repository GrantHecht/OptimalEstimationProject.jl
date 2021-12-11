function vecOuterProd!(C, a, b)
    @turbo for i ∈ eachindex(b)
        for j ∈ eachindex(a)
            C[j,i] = a[j]*b[i]
        end
    end
end
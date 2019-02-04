export GainEquality

@composite GainEquality(y, x, z, dim) begin
    @RV w = equal(x, z)
    b = uvector(dim)
    @RV y = dot(b, w)
end

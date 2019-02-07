export EqUnit

@composite EqUnit (y, x, z) begin
    @RV w = equal(x, z)
    c = [1, 0]
    @RV y = dot(c, w)
end

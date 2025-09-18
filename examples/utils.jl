using Plots

ENV["GKSwstype"]="nul"

function mismip_plus_bed(x,y)
    xbar = 300000.0
    b0 = -150.0; b2 = -728.8; b4 = 343.91; b6 = -50.75;
    wc = 24000.0; fc = 4000.0; dc = 500.0;
    bx(x)= b0 + b2 * (x / xbar) ^ 2 + b4 * (x / xbar) ^ 4 + b6 * (x / xbar) ^ 6;
    by(y)= dc *( (1 + exp(-2(y - wc) / fc)) ^ (-1) + (1 + exp(2(y + wc) / fc)) ^ (-1) );
    b = max(bx(x) + by(y), -720.0);
    return b;
end

function plot_field(data, name="test.png")
    # Now we can plot
    Plots.heatmap(data);
    plot!(size = (600, 400), fmt=:png)
    savefig(name)
end
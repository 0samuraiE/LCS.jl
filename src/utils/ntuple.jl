macro ntuple(value, n)
    make = quote
        ntuple(_ -> $value, $n)
    end
    esc(make)
end

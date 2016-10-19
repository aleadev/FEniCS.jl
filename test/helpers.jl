macro inferred_return_type(expr)
  types = Base.gen_call_with_extracted_types(Base.return_types, expr)
  quote
    @assert length($types)==1
    $types[1]
  end
end

#macro inferred_isa(expr)

"""
Impertinent copy of `Base.inferred` with type equality replaced with `isa`
"""
macro inferred_isa(ex)
    if Meta.isexpr(ex, :ref)
        ex = Expr(:call, :getindex, ex.args...)
    end
    Meta.isexpr(ex, :call)|| error("@inferred_isa requires a call expression")

    Base.remove_linenums!(quote
        $(if any(a->(Meta.isexpr(a, :kw) || Meta.isexpr(a, :parameters)), ex.args)
            # Has keywords
            args = gensym()
            kwargs = gensym()
            quote
                $(esc(args)), $(esc(kwargs)), result = $(esc(Expr(:call, _args_and_call, ex.args[2:end]..., ex.args[1])))
                inftypes = $(Base.gen_call_with_extracted_types(Base.return_types, :($(ex.args[1])($(args)...; $(kwargs)...))))
            end
        else
            # No keywords
            quote
                args = ($([esc(ex.args[i]) for i = 2:length(ex.args)]...),)
                result = $(esc(ex.args[1]))(args...)
                inftypes = Base.return_types($(esc(ex.args[1])), Base.typesof(args...))
            end
        end)
        @assert length(inftypes) == 1
        rettype = isa(result, Type) ? Type{result} : typeof(result)
        rettype <: inftypes[1] || error("return type $rettype is not a subtype of the inferred return type $(inftypes[1])")
        result
    end)
end

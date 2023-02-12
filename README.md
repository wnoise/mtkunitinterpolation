This is a response to https://discourse.julialang.org/t/using-interpolation-with-unit-verification-in-modelingtoolkit/93837 and https://julialang.slack.com/archives/CN04R7WKE/p1675271453652229

This is a demonstration of one way to combine:
 1. Units
 2. ModelingToolKit
 3. DataInterpolations

There are two obvious ways to combine units with the SciML ecosystem.
The first is to use `Unitful` quantities directly.  This basicly only
works in extremely simple scenarios.  Many solvers fail to generate
correctly typed intermediate arrays, especially when the types aren't
uniform.  It pretty much completely fails for MTK.

The second way is to use the annotation support of MTK to generate
systems that have already been unit-checked, but are internally
strictly unitless.  Therefore, all values that interact with the system
internally need to be unitless.  This means that any functions need to
return the right values when taking unitless numbers.  Conversely, your
functions also need to communicate their types to MTK.  If you don't
specially override `get_unit` for your function, it does something
semi-smart and just calls your function with unit arguments of the right
type, and looks at the units of the return value.

Unfortunately, while `DataInterpolations` works just fine with `Unitful`
quantities, it will not automatically handle both unitless and Unitful
quantities in the same table.  These are unit errors that it should
generally very much catch!

The obvious options are to:
 1. Use another interpolation package
 2. Wrap it in a function, defining `get_unit` for that function(-type).
 3. Wrap it in another structure, defining evaluation and `get_unit`.

Although I haven't done a survey I would expect most interpolation
packages to either work with units the same as `DataInterpolations`, or
for them to just fail.

The other two options are minor variations, but 3. is slightly less
boilerplate, as you only need to define things once, rather than for
each usage.  (And some use-cases get a bit complex: I don't know
how to automatically define `get_unit` for run-time generated functions,
though registering appears doable, but not obvious.)

There is still a bit of freedom of whether the wrapped table should have
the Unitful quantities itself, or just the numeric values, with the units
separate.  I lean towards the latter, and do that here.

In this repo are four example files:
 1. ./post.jl -- The original sample code, demonstrating fail on solving.
 2. ./post-fail-on-validation.jl -- demonstrating fail on validation.
 3. ./wrapped-function.jl -- failed attempt at wrapping with functions
 4. ./wrapped-struct.jl -- demonstrating a wrapping type.

Numbers 3 and 4 have modified `T` to be a proper @variable everywhere.
The struct wrapper is somewhat awkward in that `@register_symbolic` doesn't
appear to work for callable structs, so some magic is required.

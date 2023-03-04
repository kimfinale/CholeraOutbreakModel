### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 851c03a4-e7a4-11ea-1652-d59b7a6599f0
# setting up an empty package environment
begin
	import Pkg
	Pkg.activate(mktempdir())
	Pkg.Registry.update()
end

# ╔═╡ d6ee91ea-e750-11ea-1260-31ebf3ec6a9b
# add (ie install) a package to our environment
begin
	Pkg.add("Compose")
	# call `using` so that we can use it in our code
	using Compose
end

# ╔═╡ 5acd58e0-e856-11ea-2d3d-8329889fe16f
begin
	Pkg.add("PlutoUI")
	using PlutoUI
end

# ╔═╡ 07e3ea15-7151-4056-a85d-6574f003a706
using SymPy

# ╔═╡ 87adc2ad-2c76-4ba3-b122-97f3bb1271c6
md"## Equations for the SEIARW model

$$\frac{\textrm{d}S}{\textrm{d}t} =-\lambda \cdot S\\$$

"

# ╔═╡ 5b1e6995-37b4-4629-bac1-587a6dd2816a
md"$$\frac{\textrm{d}E}{\textrm{d}t} = \lambda S - \epsilon E\\$$"

# ╔═╡ 083e15dc-b1dc-4fec-b896-cedb63732458
md"$$\frac{\textrm{d}I}{\textrm{d}t} = \epsilon (1-f) E - \gamma I\\$$"

# ╔═╡ 6a30258c-6025-45af-808b-589d1e3a7e74
md"$$\frac{\textrm{d}A}{\textrm{d}t} = \epsilon f E - \gamma A\\$$"

# ╔═╡ a55cd5c4-c7d8-41a8-9006-0eac73b182b0
md"$$\frac{\textrm{d}R}{\textrm{d}t} =\gamma (A + I) \\$$"

# ╔═╡ 9d73d89b-67f6-4caa-8a07-47979b06256e
md"$$\frac{\textrm{d}W}{\textrm{d}t} =\kappa \cdot \left(I + rA\right)  - \xi W\\$$"

# ╔═╡ dccc7df6-2074-4bc6-9595-ff92cac16b66
md"$${\lambda} = \beta\left(I + rA\right)  + \beta_W \frac{W}{K+W} \\$$"

# ╔═╡ f2bbea46-a6f7-4136-8804-f53b991bbd8a
md"## Computing $R_0$ for the SEIARCW model

This is just to record the calculation required for the basic reproductive number ($\mathcal{R}_0$) for the SIRCW model. 
We used the next generation matrix method of van den Driessche and Watmough 
to derive the expression for $\mathcal{R}_0$ as a function of the model parameters.

We follow the approach by van den Driessche and Watmough (2002) Math Biosci. Let $\mathcal{F}_i$ be the rate of appearance of new infections in compartment $i$, $\mathcal{V}_i^{+}$ be the rate of transfer of individuals into compartment $i$ by all other means, and $\mathcal{V}_i^{-}$ be the rate of transfer of individuals out of compartment $i$."


# ╔═╡ 8099eff9-a139-42f7-8b46-22500a0c6586
md"The disease transmission model can be defiend as $$\dot{x} = f(x) = \mathcal{F} - \mathcal{V}$$ where $\mathcal{V} = \mathcal{V}^{-} -\mathcal{V}^{+}$"

# ╔═╡ 54412a6b-616a-470d-8dbd-d73ea75e7246
md"
$$\mathrm{\mathbf{X}} =\begin{bmatrix}
    \dot{I}\\
    \dot{A}\\
	\dot{W}\\
    \end{bmatrix} = \mathbf{F} - \mathbf{V}$$
"

# ╔═╡ 81ddbabb-50d3-4ee8-a3c0-18e618e9b648
md"
$$\mathbf{F} = \frac{\partial \mathcal{F}_i}{\partial \mathcal{X}_j} = \begin{bmatrix}
    \beta& r \beta& \frac {\beta_W K N}{(K+W)^2}\\
    0& 0& 0\\
    0& 0& 0\\
	\end{bmatrix}$$ 
"

# ╔═╡ 60ebdec8-b34b-42f9-8c4f-8ff1dffb326e
md"
$$\frac {\mathrm{d}\beta_W \frac{N W}{K+W}}{dW} = \frac {\beta_W K N}{(K+W)^2}= \frac {\beta_W N}{K} \\$$ for $W=0$ at the disease-free equilibrium"

# ╔═╡ 0387bef8-964a-4f27-b04e-70327c626ae9
md"
$$\mathbf{V} = \frac{\partial \mathcal{V}_i}{\partial \mathcal{X}_j}  = \begin{bmatrix}
    \gamma &0 &  0\\
    0& \gamma& 0\\
    -\kappa& -r \kappa& \xi\\
	\end{bmatrix}$$
"

# ╔═╡ 0cc57109-30c7-4811-a12c-8e1df2352da4
@vars β b γ ξ κ r N W K

# ╔═╡ 92f03229-b1ad-4452-a623-244ec9f7e482
md"b for the water-borne transmissibiliity $\beta_W$ as $\beta w$ is written as those they are two letters"

# ╔═╡ 92ea5729-72d5-43ff-b121-1bacd5313d48
F = [β r*β b*N/K; 0 0 0;  0 0 0]

# ╔═╡ c66d69b6-043d-4640-b0ad-76f326907315
V = [γ 0 0; 0 γ 0; -κ -r*κ ξ]

# ╔═╡ a65b5d32-c314-4c00-9a05-6143b322a1bb
FVinv = F*(V^(-1))

# ╔═╡ f9ac115a-b578-4fc5-aff9-acc32cb9ffd2
# eigval = FVinv.eigenvals()
eigvals = SymPy.eigvals(FVinv)

# ╔═╡ c8013c21-2b21-4eee-88b0-9871d64ace14
simplify(eigvals[1])

# ╔═╡ 6bf28899-5d30-4d33-964e-a579103a717e
simplify(eigvals[2])

# ╔═╡ e0a00754-9c76-430f-8286-35d3589e879b
md"## Vaccination rate 

For a vaccine campaign with coverage $\rho$ and a duration of $\Delta t$, vaccination rate $r$ per captita per unit time is as follows:

$\frac{-log(1-\rho)}{t}.$

This is because 

$N_t = N_0 e^{-r\Delta t} \implies \frac{N_t}{N_0} = 1 - \rho = e^{-r\Delta t},$

where $N_t, N_0$ indicates the initial population size and the size of the population who is not vaccinated after the vaccination campaign, respectively.

We can also the discrete version in which the unit time is $0 < \tau =< 1$ as follows

$N_t = N_0 (1 -r)^{\Delta t / \tau } \implies \frac{N_t}{N_0} = 1-\rho = (1-r)^{ \Delta t / \tau } = \mathrm{log} (1-\rho) = (\Delta t / \tau ) \mathrm{log}(1-r)$

$r = 1 - e^ { \frac {\mathrm{log}(1-\rho)}{\Delta t / \tau } }$
"

# ╔═╡ fafae38e-e852-11ea-1208-732b4744e4c2
md"_homework 0, version 2_"

# ╔═╡ 7308bc54-e6cd-11ea-0eab-83f7535edf25
# edit the code below to set your name and kerberos ID (i.e. email without @mit.edu)

student = (name = "Jazzy Doe", kerberos_id = "jazz")

# press the ▶ button in the bottom right of this cell to run your edits
# or use Shift+Enter

# you might need to wait until all other cells in this notebook have completed running. 
# scroll down the page to see what's up

# ╔═╡ cdff6730-e785-11ea-2546-4969521b33a7
md"""

Submission by: **_$(student.name)_** ($(student.kerberos_id)@mit.edu)
"""

# ╔═╡ a735187f-f06b-4f64-b3b9-5a96fdaca254


# ╔═╡ a2181260-e6cd-11ea-2a69-8d9d31d1ef0e
md"""
# Homework 0: Getting up and running

First of all, **_welcome to the course!_** We are excited to teach you about real world applications of scientific computing, using the same tools that we work with ourselves.

Before we start next week, we'd like everyone to **submit this zeroth homework assignment**. It will not affect your grade, but it will help us get everything running smoothly when the course starts. If you're stuck or don't have much time, just fill in your name and ID and submit 🙂
"""

# ╔═╡ 094e39c8-e6ce-11ea-131b-07c4a1199edf


# ╔═╡ 31a8fbf8-e6ce-11ea-2c66-4b4d02b41995


# ╔═╡ 339c2d5c-e6ce-11ea-32f9-714b3628909c
md"## Exercise 1 - _Square root by Newton's method_

Computing the square of a number is easy -- you just multiply it with itself.

But how does one compute the square root of a number?

##### Algorithm:

Given: $x$

Output: $\sqrt{x}$

1. Take a guess `a`
1. Divide `x` by `a`
1. Set a = the average of `x/a` and `a`. (The square root must be between these two numbers. Why?)
1. Repeat until `x/a` is roughly equal to `a`. Return `a` as the square root.

In general, you will never get to the point where `x/a` is _exactly_ equal to `a`. So if our algorithm keeps going until `x/a == a`, then it will get stuck.

So instead, the algorithm takes a parameter `error_margin`, which is used to decide when `x/a` and `a` are close enough to halt.
"

# ╔═╡ 56866718-e6ce-11ea-0804-d108af4e5653
md"### Exercise 1.1

Step 3 in the algorithm sets the new guess to be the average of `x/a` and the old guess `a`.

This is because the square root must be between the numbers `x/a` and `a`. Why?
"

# ╔═╡ bccf0e88-e754-11ea-3ab8-0170c2d44628
ex_1_1 = md"""
your answer here
""" 

# you might need to wait until all other cells in this notebook have completed running. 
# scroll down the page to see what's up

# ╔═╡ e7abd366-e7a6-11ea-30d7-1b6194614d0a
if !(@isdefined ex_1_1)
	md"""Do not change the name of the variable - write you answer as `ex_1_1 = "..."`"""
end

# ╔═╡ d62f223c-e754-11ea-2470-e72a605a9d7e
md"### Exercise 1.2

Write a function newton_sqrt(x) which implements the above algorithm."

# ╔═╡ 4896bf0c-e754-11ea-19dc-1380bb356ab6
function newton_sqrt(x, error_margin=0.01, a=x / 2) # a=x/2 is the default value of `a`
	return x # this is wrong, write your code here!
end

# ╔═╡ 7a01a508-e78a-11ea-11da-999d38785348
newton_sqrt(2)

# ╔═╡ 682db9f8-e7b1-11ea-3949-6b683ca8b47b
let
	result = newton_sqrt(2, 0.01)
	if !(result isa Number)
		md"""
!!! warning "Not a number"
    `newton_sqrt` did not return a number. Did you forget to write `return`?
		"""
	elseif abs(result - sqrt(2)) < 0.01
		md"""
!!! correct
    Well done!
		"""
	else
		md"""
!!! warning "Incorrect"
    Keep working on it!
		"""
	end
end

# ╔═╡ 088cc652-e7a8-11ea-0ca7-f744f6f3afdd
md"""
!!! hint
    `abs(r - s)` is the distance between `r` and `s`
"""

# ╔═╡ c18dce7a-e7a7-11ea-0a1a-f944d46754e5
md"""
!!! hint
    If you're stuck, feel free to cheat, this is homework 0 after all 🙃

    Julia has a function called `sqrt`
"""

# ╔═╡ 5e24d95c-e6ce-11ea-24be-bb19e1e14657
md"## Exercise 2 - _Sierpinksi's triangle_

Sierpinski's triangle is defined _recursively_:

- Sierpinski's triangle of complexity N is a figure in the form of a triangle which is made of 3 triangular figures which are themselves Sierpinski's triangles of complexity N-1.

- A Sierpinski's triangle of complexity 0 is a simple solid equilateral triangle
"

# ╔═╡ 6b8883f6-e7b3-11ea-155e-6f62117e123b
md"To draw Sierpinski's triangle, we are going to use an external package, [_Compose.jl_](https://giovineitalia.github.io/Compose.jl/latest/tutorial). Let's set up a package environment and add the package.

A package contains a coherent set of functionality that you can often use a black box according to its specification. There are [lots of Julia packages](https://juliahub.com/ui/Home).
"

# ╔═╡ dbc4da6a-e7b4-11ea-3b70-6f2abfcab992
md"Just like the definition above, our `sierpinksi` function is _recursive_: it calls itself."

# ╔═╡ 02b9c9d6-e752-11ea-0f32-91b7b6481684
complexity = 3

# ╔═╡ 1eb79812-e7b5-11ea-1c10-63b24803dd8a
if complexity == 3 
	md"""
Try changing the value of **`complexity` to `5`** in the cell above. 

Hit `Shift+Enter` to affect the change.
	"""
else
	md"""
**Great!** As you can see, all the cells in this notebook are linked together by the variables they define and use. Just like a spreadsheet!
	"""
end

# ╔═╡ d7e8202c-e7b5-11ea-30d3-adcd6867d5f5
md"### Exercise 2.1

As you can see, the total area covered by triangles is lower when the complexity is higher."

# ╔═╡ f22222b4-e7b5-11ea-0ea0-8fa368d2a014
md"""
Can you write a function that computes the _area of `sierpinski(n)`_, as a fraction of the area of `sierpinski(0)`?

So:
```
area_sierpinski(0) = 1.0
area_sierpinski(1) = 0.??
...
```
"""

# ╔═╡ ca8d2f72-e7b6-11ea-1893-f1e6d0a20dc7
function area_sierpinski(n)
	return 1.0
end

# ╔═╡ 71c78614-e7bc-11ea-0959-c7a91a10d481
if area_sierpinski(0) == 1.0 && area_sierpinski(1) == 3 / 4
	md"""
!!! correct
    Well done!
	"""
else
	md"""
!!! warning "Incorrect"
    Keep working on it!
	"""
end

# ╔═╡ c21096c0-e856-11ea-3dc5-a5b0cbf29335
md"**Let's try it out below:**"

# ╔═╡ 52533e00-e856-11ea-08a7-25e556fb1127
md"Complexity = $(@bind n Slider(0:6, show_value=true))"

# ╔═╡ c1ecad86-e7bc-11ea-1201-23ee380181a1
md"""
!!! hint
    Can you write `area_sierpinksi(n)` as a function of `area_sierpinski(n-1)`?
"""

# ╔═╡ c9bf4288-e6ce-11ea-0e13-a36b5e685998


# ╔═╡ a60a492a-e7bc-11ea-0f0b-75d81ce46a01
md"That's it for now, see you next week!"

# ╔═╡ b3c7a050-e855-11ea-3a22-3f514da746a4
if student.kerberos_id === "jazz"
	md"""
!!! danger "Oops!"
    **Before you submit**, remember to fill in your name and kerberos ID at the top of this notebook!
	"""
end

# ╔═╡ d3625d20-e6ce-11ea-394a-53208540d626


# ╔═╡ dfdeab34-e751-11ea-0f90-2fa9bbdccb1e
triangle() = compose(context(), polygon([(1, 1), (0, 1), (1 / 2, 0)]))

# ╔═╡ b923d394-e750-11ea-1971-595e09ab35b5
# It does not matter which order you define the building blocks (functions) of the
# program in. The best way to organize code is the one that promotes understanding.

function place_in_3_corners(t)
	# Uses the Compose library to place 3 copies of t
	# in the 3 corners of a triangle.
	# treat this function as a black box,
	# or learn how it works from the Compose documentation here https://giovineitalia.github.io/Compose.jl/latest/tutorial/#Compose-is-declarative-1
	compose(context(),
			(context(1 / 4,   0, 1 / 2, 1 / 2), t),
			(context(0, 1 / 2, 1 / 2, 1 / 2), t),
			(context(1 / 2, 1 / 2, 1 / 2, 1 / 2), t))
end

# ╔═╡ e2848b9a-e703-11ea-24f9-b9131434a84b
function sierpinski(n)
	if n == 0
		triangle()
	else
		t = sierpinski(n - 1) # recursively construct a smaller sierpinski's triangle
		place_in_3_corners(t) # place it in the 3 corners of a triangle
	end
end

# ╔═╡ 9664ac52-e750-11ea-171c-e7d57741a68c
sierpinski(complexity)

# ╔═╡ df0a4068-e7b2-11ea-2475-81b237d492b3
sierpinski.(0:6)

# ╔═╡ 147ed7b0-e856-11ea-0d0e-7ff0d527e352
md"""

Sierpinski's triangle of complexity $(n)

 $(sierpinski(n))

has area **$(area_sierpinski(n))**

"""

# ╔═╡ Cell order:
# ╠═87adc2ad-2c76-4ba3-b122-97f3bb1271c6
# ╠═5b1e6995-37b4-4629-bac1-587a6dd2816a
# ╠═083e15dc-b1dc-4fec-b896-cedb63732458
# ╠═6a30258c-6025-45af-808b-589d1e3a7e74
# ╠═a55cd5c4-c7d8-41a8-9006-0eac73b182b0
# ╠═9d73d89b-67f6-4caa-8a07-47979b06256e
# ╠═dccc7df6-2074-4bc6-9595-ff92cac16b66
# ╠═f2bbea46-a6f7-4136-8804-f53b991bbd8a
# ╠═8099eff9-a139-42f7-8b46-22500a0c6586
# ╠═54412a6b-616a-470d-8dbd-d73ea75e7246
# ╠═81ddbabb-50d3-4ee8-a3c0-18e618e9b648
# ╠═60ebdec8-b34b-42f9-8c4f-8ff1dffb326e
# ╠═0387bef8-964a-4f27-b04e-70327c626ae9
# ╠═07e3ea15-7151-4056-a85d-6574f003a706
# ╠═0cc57109-30c7-4811-a12c-8e1df2352da4
# ╠═92f03229-b1ad-4452-a623-244ec9f7e482
# ╠═92ea5729-72d5-43ff-b121-1bacd5313d48
# ╠═c66d69b6-043d-4640-b0ad-76f326907315
# ╠═a65b5d32-c314-4c00-9a05-6143b322a1bb
# ╠═f9ac115a-b578-4fc5-aff9-acc32cb9ffd2
# ╠═c8013c21-2b21-4eee-88b0-9871d64ace14
# ╠═6bf28899-5d30-4d33-964e-a579103a717e
# ╠═e0a00754-9c76-430f-8286-35d3589e879b
# ╟─fafae38e-e852-11ea-1208-732b4744e4c2
# ╟─cdff6730-e785-11ea-2546-4969521b33a7
# ╠═7308bc54-e6cd-11ea-0eab-83f7535edf25
# ╠═a735187f-f06b-4f64-b3b9-5a96fdaca254
# ╟─a2181260-e6cd-11ea-2a69-8d9d31d1ef0e
# ╟─094e39c8-e6ce-11ea-131b-07c4a1199edf
# ╟─31a8fbf8-e6ce-11ea-2c66-4b4d02b41995
# ╟─339c2d5c-e6ce-11ea-32f9-714b3628909c
# ╟─56866718-e6ce-11ea-0804-d108af4e5653
# ╠═bccf0e88-e754-11ea-3ab8-0170c2d44628
# ╟─e7abd366-e7a6-11ea-30d7-1b6194614d0a
# ╟─d62f223c-e754-11ea-2470-e72a605a9d7e
# ╠═4896bf0c-e754-11ea-19dc-1380bb356ab6
# ╠═7a01a508-e78a-11ea-11da-999d38785348
# ╟─682db9f8-e7b1-11ea-3949-6b683ca8b47b
# ╟─088cc652-e7a8-11ea-0ca7-f744f6f3afdd
# ╟─c18dce7a-e7a7-11ea-0a1a-f944d46754e5
# ╟─5e24d95c-e6ce-11ea-24be-bb19e1e14657
# ╟─6b8883f6-e7b3-11ea-155e-6f62117e123b
# ╠═851c03a4-e7a4-11ea-1652-d59b7a6599f0
# ╠═d6ee91ea-e750-11ea-1260-31ebf3ec6a9b
# ╠═5acd58e0-e856-11ea-2d3d-8329889fe16f
# ╟─dbc4da6a-e7b4-11ea-3b70-6f2abfcab992
# ╠═e2848b9a-e703-11ea-24f9-b9131434a84b
# ╠═9664ac52-e750-11ea-171c-e7d57741a68c
# ╠═02b9c9d6-e752-11ea-0f32-91b7b6481684
# ╟─1eb79812-e7b5-11ea-1c10-63b24803dd8a
# ╟─d7e8202c-e7b5-11ea-30d3-adcd6867d5f5
# ╠═df0a4068-e7b2-11ea-2475-81b237d492b3
# ╟─f22222b4-e7b5-11ea-0ea0-8fa368d2a014
# ╠═ca8d2f72-e7b6-11ea-1893-f1e6d0a20dc7
# ╟─71c78614-e7bc-11ea-0959-c7a91a10d481
# ╟─c21096c0-e856-11ea-3dc5-a5b0cbf29335
# ╟─52533e00-e856-11ea-08a7-25e556fb1127
# ╟─147ed7b0-e856-11ea-0d0e-7ff0d527e352
# ╟─c1ecad86-e7bc-11ea-1201-23ee380181a1
# ╟─c9bf4288-e6ce-11ea-0e13-a36b5e685998
# ╟─a60a492a-e7bc-11ea-0f0b-75d81ce46a01
# ╟─b3c7a050-e855-11ea-3a22-3f514da746a4
# ╟─d3625d20-e6ce-11ea-394a-53208540d626
# ╟─dfdeab34-e751-11ea-0f90-2fa9bbdccb1e
# ╟─b923d394-e750-11ea-1971-595e09ab35b5

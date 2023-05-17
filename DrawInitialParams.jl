### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# ╔═╡ c15809e0-e919-11ed-2c6b-e581c09e359e
using Distributions ,Random ,DelimitedFiles ,BenchmarkTools

# ╔═╡ 8525a0df-90e7-42a2-b0c5-170a04544345
using Base.Threads

# ╔═╡ f4b2a2bd-bc2c-43e5-84f4-a334ee23db15
md"""
## Yadira Gaibor
## Course 8
## For 18.377 Final Project
"""

# ╔═╡ 532d4d22-595a-41f5-8c1c-b6b9166ff397
md"""
# Parallelized drawing initial parameters
"""

# ╔═╡ 1d306175-f384-45c8-bc0e-c978ac079117
begin
	#Adapted from https://python4mpia.github.io/fitting_data/MC-sampling-from-Salpeter.html
	# Draw random samples from Kroupa IMF.
	# N     ... number of samples.
	# alpha ... power-law index.
	# M_min ... lower bound of mass interval.
	# M_max ... upper bound of mass interval.
	
	function sampleFromIMF(N, alpha, M_min, M_max)
	    # Convert limits from M to logM.
	    log_M_Min = log(M_min)
	    log_M_Max = log(M_max)
	    N=1
	
	    # Since IMF decays, maximum likelihood occurs at M_min
	    maxlik = float(M_min)^ (-alpha)
	
	    # Prepare array for output masses.
	    Masses = Float64[]
	    # Fill in array.
	    while length(Masses) < N
	        # Draw candidate from logM interval.
	        logM = rand(Uniform(log_M_Min, log_M_Max))
	        M = round(exp(logM), digits=2)
	        # Compute likelihood of candidate from IMF.
	        likelihood = float(M)^(-alpha)
	        # Accept randomly.
	        u = rand(Uniform(0.0, maxlik))
	        if u < likelihood
	            push!(Masses, M)
	        end
	    end
	    return Masses[1]
	end
	
	
	
	#Draw sample of M_prog, M_star, P_prog.
	
	function draw_init_params()
	    #Draw initial period function, in days
	    mu, sigma = 4.8,2.3 # mean and standard deviation, from Kroupa et al. A&A 529, A92 (2011), and Duquennoy & Mayor (1991)
	
		s=rand(LogNormal(mu/log(10), sigma/log(10)), 1)[1] # convert to base 10
	  
	    #in days
	    minP=0.1
	    maxP=10^8
	    
	#     print('before',s)
	    while (s <minP) || (s >maxP) #make sure the period is within specified range
	        s = rand(LogNormal(mu/log(10), sigma/log(10)), 1)[1]
		end
	    
	    P_prog= s /365 #to years
	    
	    # Draw samples Mass in solar
	    M_prog = sampleFromIMF(1, 2.35, 2.0, 150.0)[1] #min:2.0, max:150
		
	    #Draw mass ratio
	    q=round(rand(Uniform(0.1, 1), 1)[1], digits=3) #Defined as M2/M1
	    
	    M_star= q * M_prog #companion mass, can read these values in from EB catalog companions
	
	    return [M_prog, q, M_star, P_prog]
	end
	
	
	
end

# ╔═╡ bac5b25e-dfe9-454e-80c3-9ecac1b3d77b
# md"""
# ### Non-parallel version
# """

# ╔═╡ 60dca83e-fd0f-4106-a25f-144448e01ba0
function non_parallel_version(Samples)
    Init_parNP = []
    for i in 1:Samples
        calc_par = draw_init_params()
        push!(Init_parNP, calc_par)
    end
    return Init_parNP
end

# ╔═╡ 2d76decb-a3bc-4edb-b717-5134cb3e52b1
begin
	# Benchmark the non-parallelized version
	Samples=[10,100,1000,10000,100000,250000]
	
	for i in 1:length(Samples)
		@time non_parallel_results = non_parallel_version(Samples[i])
	end
	
end

# ╔═╡ f19c23a5-2bc9-4fda-ad66-04c3980e1f6e
# md"""
# ### Parallel version
# """

# ╔═╡ 366662c6-3b98-40e0-9d5b-b3229ca24bde
function parallel_version(Samples)
    Init_parP = []
    @threads for i in 1:Samples
        calc_par = draw_init_params()
        push!(Init_parP, calc_par)
    end

    return Init_parP
end

# ╔═╡ f9e98f19-ec2d-444d-8d10-bb06cac621a1
# Benchmark the parallelized version
begin
	# Benchmark the non-parallelized version
	Samples2=[10,100,1000,10000,100000]#,150000]
	for i in 1:length(Samples2)
		@time parallel_results = parallel_version(Samples2[i])
	end
	
end

# ╔═╡ 74c9959e-79de-4e72-84a7-3c6de6c0e2b7
md"""
#### Save initial parameters to file
"""

# ╔═╡ be53245f-7a4d-4804-9dbc-2e5cb2ab3cad
# begin
# 	#Save samples to txt file
# 	datafile = open("Initial_binariesNP.txt", "w")
	
# 	for i in 1:length(non_parallel_results)
# 	    write(datafile, join(string.(non_parallel_results[i]), "\t") * "\n")
# 	end
# 	close(datafile)
# end

# ╔═╡ 1443d3d8-a5be-418c-bac6-563c5761e210
begin
	#Save samples to txt file
	parallel_results_save = parallel_version(1000)
	
	datafile = open("Initial_binariesNP.txt", "w")
	
	for i in 1:length(parallel_results_save)
	    write(datafile, join(string.(parallel_results_save[i]), "\t") * "\n")
	end
	close(datafile)
end

# ╔═╡ fbf1d0e5-6ab8-461c-88ac-99a19ed576e2

# begin
# 	#Save samples to txt file
# 	datafileP = open("Initial_binariesPar.txt", "w")
	
# 	for i in 1:length(parallel_results)
# 	    write(datafileP, join(string.(parallel_results[i]), "\t") * "\n")
# 	end
# 	close(datafileP)
# end


# ╔═╡ b1f4918c-a1c2-4c1c-83c1-a777763811f7
md"""
#### Get occurrence rates and save parameters to separate files
"""

# ╔═╡ 8f8d9b10-5673-4bc7-ac59-0f1ab8eb6403
#Load in samples
Init_params = readdlm("Initial_binariesNP.txt", '\t'); #Mprog,q,Mstar,Pprog

# ╔═╡ 17863b90-ed68-49f7-bce0-566ca8f6b2ad
begin
	function CEphase(mprog, mcore, sma0, mstar)
	    alph_lam = 1 # alpha lambda free parameter
	    
	    R_roche = roche(mprog, mstar, sma0) # roche radius
	    a_f = (sma0 * mcore / mprog) * ((2 * sma0 * (mprog - mcore) / (alph_lam * R_roche * mstar)) + 1) ^(-1) # final semi-major axis
	    roche_param = R_roche / sma0
	    
	    return a_f, roche_param
	end
	
	
	
	function occurence_field(m_prog, period, M_star, jacobian_CE)
	    # prog.mass in Msolar, period in days, Companion mass, mass bins, period bins, jacobian from CE phase (must be length of M_prog)
	    M_bh = m_prog # range of M_BH [M_Sun] (or linspace?)
	    P = period
	    
	    # Heaviside step fnctn.to limit stars <5 Msun
	    H = 1
	    if M_bh < 5
	        H = 0
	    end
	            
	    q = M_star / M_bh # mass ratio. This is M2/M1
	    C = (4.825494e-2) / q # normalization constant for binary fraction [Should be same size as M_bh]
	
	    k = 0.7048153 # normalization constant so that integral=1 (see Mashian & Loeb 2017, eqn below Eqn.1)
	    IBH = k * (M_bh / 0.5)^(-2.3) # IMF for black holes. **Might have to renormalize so that stars with M>20 M_sun become BH.
	    IMF = (M_star)^(-2.3) # IMF for companion (from Kroupa 2000). This can be varied later.
	
	    # Gets occurence rate for field binaries, given M_star for a range of M_bh and periods.
	    occ_field = []
	
	    occ = (C / M_bh) * (H * IBH) / IMF
	    
	    occ_each = occ / P
	
	    fin_occ = jacobian_CE .* occ_each
	    
	    return fin_occ
	end
	
	
	
	function roche(M1, M2, a)
	    # Roche Radius
	    # M1, M2 should have the same units
	    # a is semi-major axis in AU
	    q = M1 / M2 # mass ratio
	    RL = (a * 0.49 * (q^(2/3))) / ((0.6 * (q^(2/3))) + log(1 + (q^(1/3))))
	    return RL
	end
end

# ╔═╡ acd228bd-f754-4726-81c9-10050e41115c
md"""
###Optimize occurence rates calculation, and save to files.
"""

# ╔═╡ 24e802c0-d64b-4b13-b734-1854ba26283a
begin
non_CE=[] #system info for those that dont undergo CE
CE_init=[] #Initial system info for those that undergo CE
CE_final=[] #Final system info for those that undergo CE

datafileCE = open("AfterCE_binariesJ.csv", "w")
datafileCE_init = open("InitialCE_binariesJ.csv", "w")
datafile2 = open("NonCE_binariesJ.csv", "w")

for i = 1:size(Init_params)[1]
    M_prog=Init_params[i,1]
    q = Init_params[i,2] #M2/M1
    M_star=Init_params[i,3]
    P_prog=Init_params[i,4]
    
    k=0.7048153 #Normalization constant so that integral=1 (see Mashian & Loeb 2017, eqn below Eqn.1)
    IBH= k*(M_prog/0.5)^(-2.3)

    #Mass at the end of AGB stage
    M_coreAGB=(4.36e-4 *(M_prog^(5.22)) + 6.84e-2)^(1/4) # Using Eqn.66 of Hurley et al 2000

    #AGB max radius
    b33=2.474401
    b31=7.425137^b33
    b32=9.268325
    b34=1.127018^b33

    L=(b31+b32* M_prog^(b33+1.8))/(b34+ M_prog^b33)
    

    R_max=1.125* M_prog^(-0.33) * (L^0.4 +0.383*(L^0.76)) #R in the AGB (per Hurley et al 2000), in solar radii.

    R_max=R_max*0.00465047 #Convert R_sun to AU.


    #Convert initial period in years to sma in AU
    sma_all=(P_prog^2 *(M_prog+M_star))^(1/3)
    
    #Selects systems where sma<Rmax (those undergo CE evol.)
    #For a given M_star, M_prog, and Period
    #Will go through CE if roche_radius<R_max
    sma_init_CE=0
    
    roche_radius=roche(M_prog,M_star,sma_all)

    if roche_radius<=R_max #Undergoes CE
        sma_init_CE=sma_all #just sma value for the rest of the calculations
        
        push!(CE_init,[M_prog,q,M_star,P_prog*365,sma_init_CE]) #all initial system info in days and au
		write(datafileCE_init,"$M_prog,$q,$M_star,$(P_prog*365),$sma_init_CE\n")

        
        #Undergoes CE
        roche_hold=0
        delta_sma=0

        new_sma,roche_rad=CEphase(M_prog,M_coreAGB,sma_init_CE,M_star)

        roche_hold= copy(roche_rad)
        
        delta_a=new_sma-sma_init_CE
        delta_sma= copy(delta_a)   

        M_prog_finCE=M_coreAGB #BH mass at the end of CE is just stripped core.

        #Convert sma (au) to period (years)
        period_CE=(((new_sma)^3)/(M_coreAGB+M_star))^(1/2) #in years
        period_CE=period_CE*365 #to days, still have to round to sig. figs/ give uncertainty
        
        #jacobian for CE phase
        k=0.2 #free param later
        r_L=roche_hold
        q2=M_prog/M_star #M1/M2 for this q

		alphlam=1 #free parameter later
        
        Jac_ce= ((2*(1-k))/(alphlam*r_L*k* 1/q) + 1/k)^(-1) 
        
        occ_rate=occurence_field(M_prog,P_prog*365,M_star,Jac_ce) #gets occurence rates
        
        push!(CE_final,[M_prog_finCE,q,M_star,period_CE,new_sma,occ_rate])

		write(datafileCE,"$M_prog_finCE,$q,$M_star,$period_CE,$new_sma, $occ_rate\n")
    
    else #does not undergo CE
        #jacobian
        k=0.2
        r_L=roche(M_prog,M_star,sma_all)
        q2=M_prog/M_star #M1/M2 for this q
		alphlam=1
        
        Jac=1/k*((2*(1-k)/(alphlam*r_L*k*q))+1/k)
        
        occ_rate=occurence_field(M_prog,P_prog*365,M_star,Jac) #gets occurence rates
    
        push!(non_CE,[M_prog,q,M_star,P_prog*365,sma_all,occ_rate]) #all system info in days and au
		write(datafile2,"$M_prog,$q,$M_star,$(P_prog*365),$sma_all,$occ_rate\n")
	end
end
	
close(datafileCE)
close(datafileCE_init)
close(datafile2)

end

# ╔═╡ d2b8fe2a-712d-45c5-90fc-e60d8dd2b118


# ╔═╡ deaa0777-82c0-42cf-ab75-950e8abef39b


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
DelimitedFiles = "8bb1440f-4735-579b-a4ab-409b98df4dab"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[compat]
BenchmarkTools = "~1.3.2"
Distributions = "~0.25.88"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.5"
manifest_format = "2.0"
project_hash = "03e18dea669acc4accb537f8ea4a7bef5c1ba8a3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "d9a9701b899b30332bbcb3e1679c41cce81fb0e8"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.3.2"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c6d890a52d2c4d55d326439580c3b8d0875a77d9"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.7"

[[deps.ChangesOfVariables]]
deps = ["LinearAlgebra", "Test"]
git-tree-sha1 = "f84967c4497e0e1955f9a582c232b02847c5f589"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.7"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "7a60c856b9fa189eb34f5f8a6f6b5529b7942957"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.6.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.1+0"

[[deps.DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "180538ef4e3aa02b01413055a7a9e8b6047663e1"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.88"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "fc86b4fd3eff76c3ce4f5e96e2fdfa6282722885"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.0.0"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "432b5b03176f8182bd6841fbfc42c718506a2d5f"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.15"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "6667aadd1cdee2c6cd068128b3d226ebc4fb0c67"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.9"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "0a1b7c2863e44523180fdb3146534e265a91870b"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.23"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "d321bf2de576bf25ec4d3e4360faca399afca282"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "67eae2738d63117a196f497d7db789821bce61d1"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.17"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7302075e5e06da7d000d9bfa055013e3e85578ca"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.9"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "2e47054ffe7d0a8872e977c0d09eb4b3d162ebde"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.0.2"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "6ec7ac8412e83d57e313393220879ede1740f9ee"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.8.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "ef28127915f4229c971eb43f3fc075dd3fe91880"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.2.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "45a7769a04a3cf80da1c1c7c60caf932e6f4c9f7"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.6.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "f625d686d5a88bcd2b15cd81f18f98186fdc0c9a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╠═c15809e0-e919-11ed-2c6b-e581c09e359e
# ╠═f4b2a2bd-bc2c-43e5-84f4-a334ee23db15
# ╟─532d4d22-595a-41f5-8c1c-b6b9166ff397
# ╠═8525a0df-90e7-42a2-b0c5-170a04544345
# ╠═1d306175-f384-45c8-bc0e-c978ac079117
# ╠═bac5b25e-dfe9-454e-80c3-9ecac1b3d77b
# ╠═60dca83e-fd0f-4106-a25f-144448e01ba0
# ╠═2d76decb-a3bc-4edb-b717-5134cb3e52b1
# ╠═f19c23a5-2bc9-4fda-ad66-04c3980e1f6e
# ╠═366662c6-3b98-40e0-9d5b-b3229ca24bde
# ╠═f9e98f19-ec2d-444d-8d10-bb06cac621a1
# ╠═74c9959e-79de-4e72-84a7-3c6de6c0e2b7
# ╠═be53245f-7a4d-4804-9dbc-2e5cb2ab3cad
# ╠═1443d3d8-a5be-418c-bac6-563c5761e210
# ╠═fbf1d0e5-6ab8-461c-88ac-99a19ed576e2
# ╠═b1f4918c-a1c2-4c1c-83c1-a777763811f7
# ╠═8f8d9b10-5673-4bc7-ac59-0f1ab8eb6403
# ╠═17863b90-ed68-49f7-bce0-566ca8f6b2ad
# ╠═acd228bd-f754-4726-81c9-10050e41115c
# ╠═24e802c0-d64b-4b13-b734-1854ba26283a
# ╠═d2b8fe2a-712d-45c5-90fc-e60d8dd2b118
# ╠═deaa0777-82c0-42cf-ab75-950e8abef39b
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

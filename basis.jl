using SpecialFunctions
using DelimitedFiles
using Dates

function get_bonds_v(Nx::Int, Ny::Int)
    bonds = zeros(Int, Nx*Ny*3, 4)
    for i in 0:Nx-1
        for j in 0:Ny-1
            ij = i + j * Nx
            imj = (i - 1 + Nx) % Nx + j * Nx
            ijm = i + ((j - 1 + Ny) % Ny) * Nx
            bonds[3*ij+1, :] = [6*ij+1, 6*ij+2, 6*imj+5, 6*ijm+4]
            bonds[3*ij+2, :] = [6*ij+2, 6*ij+3, 6*ij+5, 6*ijm+6]
            bonds[3*ij+3, :] = [6*ij+1, 6*ij+3, 6*ij+4, 6*imj+6]
        end
    end
    return bonds
end

function Z_gen(Nx::Int, Ny::Int)
    if Nx < 2 || Ny < 2
        throw("Nx, Ny must be at least 2")
    end
    bonds = zeros(Int, 12)
    bonds[1] = 1
    bonds[2] = 2
    bonds[3] = (Ny - 1) * Nx * 6 + 6
    bonds[4] = (1 + (Ny - 1) * Nx) * 6 + 4
    bonds[5] = 8
    bonds[6] = 9
    bonds[7] = 10
    bonds[8] = Nx * 6 + 5
    bonds[9] = Nx * 6 + 3
    bonds[10] = Nx * 6 + 1
    bonds[11] = (2 * Nx - 1) * 6 + 5
    bonds[12] = (Nx - 1) * 6 + 6
    return bonds
end

function X_gen(Nx::Int, Ny::Int)
    vert = zeros(Int, 6)
    vert[1] = 2
    vert[2] = 4
    vert[3] = 6
    vert[4] = Nx * 3 + 2
    vert[5] = Nx * 3 + 1
    vert[6] = 3
    return vert
end

function comb(a::Int, b::Int)
    return exp(logfactorial(a) - logfactorial(b) - logfactorial(a-b))
end

function get_bonds_b(Nx::Int, Ny::Int)
    bonds = Vector{Int}[]
    for j in 0:Ny-1
        for i in 0:Nx-1
            ij = 6 * (i + j * Nx)
            imj = 6 * ((i - 1 + Nx) % Nx + j * Nx)
            ipj = 6 * ((i + 1) % Nx + j * Nx)
            ijm = 6 * (i + ((j - 1 + Ny) % Ny) * Nx)
            ijp = 6 * (i + ((j + 1) % Ny) * Nx)
            imjp = 6 * ((i - 1 + Nx) % Nx + ((j + 1) % Ny) * Nx)
            ipjm = 6 * ((i + 1) % Nx + ((j - 1 + Ny) % Ny) * Nx)

            t = [ij+2, ij+3, ij+4, imj+6, imj+5, ijm+4]
            push!(bonds, ∩(setdiff(1:6*Nx*Ny, t), ij+2:6*Nx*Ny))
            t = [ij+1, ij+3, ij+5, imj+5, ijm+4, ijm+6]
            push!(bonds, ∩(setdiff(1:6*Nx*Ny, t), ij+3:6*Nx*Ny))
            t = [ij+1, ij+2, ij+4, ij+5, imj+6, ijm+6]
            push!(bonds, ∩(setdiff(1:6*Nx*Ny, t), ij+4:6*Nx*Ny))
            t = [ij+1, ij+3, imj+6, ijp+1, ijp+2, imjp+5]
            push!(bonds, ∩(setdiff(1:6*Nx*Ny, t), ij+5:6*Nx*Ny))
            t = [ij+2, ij+3, ipj+1, ipj+2, ijm+6, ipjm+4]
            push!(bonds, ∩(setdiff(1:6*Nx*Ny, t), ij+6:6*Nx*Ny))
            t = [ijp+2, ijp+3, ijp+5, ipj+1, ipj+3, ipj+4]
            push!(bonds, ∩(setdiff(1:6*Nx*Ny, t), ij+7:6*Nx*Ny))
        end
    end
    return bonds
end

function basis(bonds::Vector{Int}, bonds_dic::Vector{Vector{Int}},
        ret::Vector{Int})
    #global basis_all
    if bonds == []
        #println(ret)
        #push!(basis_all,ret)
        global basis_count += 1
        if basis_count % 1000000 == 0
            println(basis_count)
            flush(stdout)
        end
        i = length(ret) + 1
        conf_count[i] += 1
        conf = zeros(Int, 6*Nx*Ny)
        conf[ret] .= 1
        
        zc = 0
        for l in 1:12
            zc += conf[Z[l]]
        end
        conf_zc[i] += 1 - zc % 2

        zo = 0
        for l in 1:6
            zo += conf[Z[l]]
        end
        conf_zo[i] += 1 - zo % 2

        xc = 1
        for v in 1:6
            if (conf[bonds_v[X[v],1]] + conf[bonds_v[X[v],2]]
                + conf[bonds_v[X[v],3]] + conf[bonds_v[X[v],4]]) == 0
                xc = 0
                break
            end
        end
        conf_xc[i] += xc

        xo = 1
        if (conf[bonds_v[2,1]] + conf[bonds_v[2,2]]
            + conf[bonds_v[2,3]] + conf[bonds_v[2,4]] == 1
            && conf[bonds_v[3,1]] + conf[bonds_v[3,2]]
            + conf[bonds_v[3,3]] + conf[bonds_v[3,4]] == 1
            && conf[7] == 0 && conf[8] == 0 && conf[Nx*6+1] == 0
            && conf[Nx*6+2] == 0)
            xo6 = (conf[Z[1]] + conf[Z[2]] + conf[Z[3]] + conf[Z[4]]
                   + conf[Z[11]] + conf[Z[12]])
            if xo6 % 2 == 1
                conf_xo0[i] += 1
            else
                xo3 = conf[3] + conf[4] + conf[5]
                if (xo6 ÷ 2 + xo3) % 2 == 0
                    conf_xom[i] += 1
                else
                    conf_xop[i] += 1
                end
            end
        end
    else
        basis(Int[],bonds_dic,ret)
        for i in bonds
            #println(i)
            #println(vcat(ret,[i]))
            basis(∩(bonds_dic[i], bonds),bonds_dic,vcat(ret,[i]))
        end
    end
end

function main1(Nx::Int, Ny::Int)
    global bonds_v = get_bonds_v(Nx, Ny)
    global Z = Z_gen(Nx, Ny)
    global X = X_gen(Nx, Ny)
    conf_max = Nx * Ny * 6 ÷ 4 + 1
    #global basis_all = Vector{Int}[]
    global basis_count = 0
    global conf_count = zeros(Int, conf_max)
    global conf_zc = zeros(Int, conf_max)
    global conf_zo = zeros(Int, conf_max)
    global conf_xc = zeros(Int, conf_max)
    global conf_xo0 = zeros(Int, conf_max)
    global conf_xop = zeros(Int, conf_max)
    global conf_xom = zeros(Int, conf_max)
    basis(collect(1:6*Nx*Ny), get_bonds_b(Nx, Ny), Int[])
    #println(basis_all[1:20])
    #println(length(basis_all))
    println(conf_count)
    println(conf_zc)
    println(conf_zo)
    println(conf_xc)
    println(conf_xo0)
    println(conf_xop)
    println(conf_xom)
end


Nx = parse(Int,ARGS[1])
Ny = parse(Int,ARGS[2])
#conf_max = Nx * Ny * 6 ÷ 4
#coe = comb.(6*Nx*Ny,1:conf_max)/1e6

println(Dates.format(now(), "u d HH:MM:SS"))
main1(Nx, Ny)
println(Dates.format(now(), "u d HH:MM:SS"))
#=for i in 1:10
    conf_count, conf_zc, conf_zo, conf_xc, conf_xo0, conf_xop, conf_xom = main(Nx, Ny, Int(1e6))
    println(conf_count)
    println(conf_zc)
    println(conf_zo)
    println(conf_xc)
    println(conf_xo0)
    println(conf_xop)
    println(conf_xom)
    #println(conf_count.*coe)
    #println(conf_z.*coe)
    #println(conf_x.*coe)
    flush(stdout)
end=#
#println(conf_count)

using DelimitedFiles,PlotlyJS

m = readdlm("out.txt")
m = Float64.(m[2:end,:])

input = readlines("in.txt")
nc = eval(Meta.parse(input[6]))
ns = eval(Meta.parse(input[7]))

xy = reshape(eval(Meta.parse("["*input[9]*"]")),2,nc)
xys= reshape(eval(Meta.parse("["*input[10]*"]")),2,ns)
ds = eval(Meta.parse("["*input[11]*"]"))

trace1 = scatter(x=[m[1:201,5];m[20002:20202,5]],y=[m[1:201,4];m[20002:20202,4]],mode = "markers",text =  "ϵ = ".*string.([m[1:201,1];m[20002:20202,1]]).*"<br>+ ".*string.([m[1:201,2];m[20002:20202,2]]).*"y",visible  = false);

trace2 = scatter3d(x=m[:,5],y=m[:,6],z=m[:,4],mode = "markers",text =  "ϵ = ".*string.(m[:,1]).*"<br>+ ".*string.(m[:,2]).*"y".*"<br>+ ".*string.(m[:,3]).*"z",marker_sizeref = 1.5,marker_size = 1.2,opacity = 0.5,visible  = false);

trace31 = scatter(x=[xy[1,:];xy[1,1]],y=[xy[2,:];xy[2,1]],mode = "lines",visible  = true)
trace32 = scatter(x = xys[1,:],y = xys[2,:], mode = "markers",text = "d = ".*string.(ds),visible  = true)


trace = [trace31,trace32,trace1,trace2]


menu =
    attr(showactive = true,
         buttons=[
             attr(label="Prerez",
                  method="update",
                  args=[Dict(
                      "visible"=> [true,true,false,false],
                      #"scene" => attr(visible = false),
                      #"xaxis" => attr(visible=true),
                      #"yaxis" => attr(visible=true),
                      )]),
                  attr(label="Enoosni upogib",method="update", args=[Dict("visible"=> [false,false,true,false],
                                                                          #"scene" => (visible = false),
                                                                          #"xaxis" => attr(visible=true),
                                                                          #"yaxis" => attr(visible=true)
                                                                          )]),
                  attr(label="Dvoosni upogib",method="update", args=[Dict("visible"=> [false,false,false,true],
                                                                          #"scene" => attr(xaxis=attr(title="X"), yaxis=attr(title="Y"), zaxis=attr(title="Z")),
                                                                          #"xaxis" => attr(visible=false),
                                                                          #"yaxis" => attr(visible=false)
                                                                          )])])
lyout =Layout(

        #=[
        attr(showactive = true,
             buttons=[
                 attr(label="Prerez",
                      method="update",
                      args=[Dict(
                        "visible"=> [true,true,false,false],
                        #"scene" => attr(visible = false),
                        #"xaxis" => attr(visible=true),
                        #"yaxis" => attr(visible=true),
                        )]),
                 attr(label="Enoosni upogib",method="update", args=[Dict("visible"=> [false,false,true,false],
                     #"scene" => (visible = false),
                     #"xaxis" => attr(visible=true),
                     #"yaxis" => attr(visible=true)
                     )]),
            attr(label="Dvoosni upogib",method="update", args=[Dict("visible"=> [false,false,false,true],
                #"scene" => attr(xaxis=attr(title="X"), yaxis=attr(title="Y"), zaxis=attr(title="Z")),
                #"xaxis" => attr(visible=false),
                #"yaxis" => attr(visible=false)
                )])])]=#

        xaxis=attr(title="X"),
        yaxis=attr(title="Y"),
        scene=attr(xaxis=attr(title="X"),yaxis=attr(title="Y"),zaxis=attr(title="Z")),
        hovermode="closest",
        updatemenus=[menu]
        )

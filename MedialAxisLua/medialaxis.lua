label = "Medial Axis"

about = [[
Computes the medial axis of a simple convex polygon.
Select a polygon and run to generate its medial axis.
]]

ipelet = false

function run(model, num)
  if not ipelet then
    ipelet = assert(ipe.Ipelet(dllname))
  end
  model:runIpelet("Compute Medial Axis", ipelet, num)
end



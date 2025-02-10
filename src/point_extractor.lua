-- Ipelet Information
label = "Coordinate Extractor"
about = "Enumerates coordinates of selected points and polygons in Ipe"
log_file = nil

local function incorrect_input(model)
    model:warning("Please select mark (point) or path (polygon) objects!")
end

local function get_coordiantes_file()
    log_file = ipe.openFile("coordinates.txt", "w")
end

local function close_coordinates_file()
    if log_file then log_file:close() end
end

local function write_coordinates(message)
    if log_file then
        log_file:write(message .. "\n")
    else
        error("Could not open coordinates.txt for writing")
    end
end

local function check_unique_point(px, py, points)
    local key = string.format("%.2f,%.2f", px, py)
    if not points[key] then
        points[key] = true
        write_coordinates(string.format("POINT (%.2f %.2f)", px, py))
    end
end

local function check_unique_polygon_point(px, py, points, polygon_coords)
    local key = string.format("%.2f,%.2f", px, py)
    if not points[key] then
        points[key] = true
        table.insert(polygon_coords, string.format("%.2f %.2f", px, py))
    end
end

local function is_mark(obj)
    return obj:type() == "reference" and obj:symbol():sub(1, 5) == "mark/"
end

local function is_path(obj)
    return obj:type() == "path"
end

function run(model)
    local p = model:page()
    get_coordiantes_file()

    if #model:selection() == 0 then
        model:warning("No objects selected!")
        return
    end

    for _, obj, sel, layer in p:objects() do
        if sel and is_mark(obj) then
            -- Extract and log unique point coordinates
            local point = obj:matrix() * obj:position()
            check_unique_point(point.x, point.y, {})
        elseif sel and is_path(obj) then
            -- Extract and log polygon coordinates
            local matrix = obj:matrix()
            local shape = obj:shape()

            local unique_points = {}  -- Track unique points
            local polygon_coords = {} -- Store polygon coordinates
            local first_point = nil   -- Store first point for closure

            for _, subpath in ipairs(shape) do
                if subpath.type == "closed" then
                    for _, point in ipairs(subpath) do
                        local transformed_point = matrix * point
                        check_unique_polygon_point(transformed_point.x, transformed_point.y, unique_points,
                            polygon_coords)

                        -- Store first point for closure
                        if not first_point then
                            first_point = string.format("%.2f %.2f", transformed_point.x, transformed_point.y)
                        end
                    end
                elseif subpath.type == "curve" then
                    for _, segment in ipairs(subpath) do
                        for _, point in ipairs(segment) do
                            local transformed_point = matrix * point
                            check_unique_polygon_point(transformed_point.x, transformed_point.y, unique_points,
                                polygon_coords)

                            -- Store first point for closure
                            if not first_point then
                                first_point = string.format("%.2f %.2f", transformed_point.x, transformed_point.y)
                            end
                        end
                    end
                end
            end

            -- Ensure the polygon is closed by adding the first point again
            if first_point then
                table.insert(polygon_coords, first_point)
            end

            -- Write the WKT representation of the polygon
            write_coordinates("POLYGON ((" .. table.concat(polygon_coords, ", ") .. "))")
        elseif sel then
            incorrect_input(model)
        end
    end

    model:warning("Coordinates saved to coordinates.txt!")
    close_coordinates_file()
end

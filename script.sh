#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# =========================
# Set de opciones a correr
# Puedes sobreescribir cada variable por ENV al ejecutar el script.
# Ejemplo:
# FAIRNESS_SET="NONE,PEA,SA" INSTANCE_FROM=1 INSTANCE_TO=5 ./script.sh
# =========================

DEFAULT_INST_FOLDER="inst/inst2020"
if [[ ! -d "$DEFAULT_INST_FOLDER" && -d "../inst/inst2020" ]]; then
  DEFAULT_INST_FOLDER="../inst/inst2020"
fi
SCRIPT_NAME="${SCRIPT_NAME:-script.sh}"
INST_FOLDER="${INST_FOLDER:-$DEFAULT_INST_FOLDER}"
INSTANCE_FROM="${INSTANCE_FROM:-1}"
INSTANCE_TO="${INSTANCE_TO:-10}"

# Medidas disponibles: NONE, PEA, SA, LEXMMFPEA, LEXMMFSA
# Este script es el motor genérico y requiere FAIRNESS_SET explícito.
FAIRNESS_SET="${FAIRNESS_SET:-}"
if [[ -z "$FAIRNESS_SET" ]]; then
  echo "ERROR: define FAIRNESS_SET (ejemplo: FAIRNESS_SET=PEA)" >&2
  exit 1
fi

# Formato: "NBstage:childs:periods;NBstage:childs:periods;..."
TREE_SET="${TREE_SET:-6:2:4,6:6:4,8:2:3,8:6:3}"

# Sets de parametros numéricos (coma separada)
J_SET="${J_SET:-5,10}"
THETA_SET="${THETA_SET:-0.2,0.6}"
AVG_D_SET="${AVG_D_SET:-100.0}"
DEV_D_SET="${DEV_D_SET:-10.0,20.0}"

OUT_CSV="${OUT_CSV:-fairness_config_set_report.csv}"
OUT_CSV_HOUSE="${OUT_CSV_HOUSE:-fairness_config_set_report_by_house.csv}"
LOCK_OUTPUT="${LOCK_OUTPUT:-1}"

if [[ "$LOCK_OUTPUT" == "1" ]]; then
  if ! command -v flock >/dev/null 2>&1; then
    echo "ERROR: 'flock' no está disponible y LOCK_OUTPUT=1. Instala flock o usa LOCK_OUTPUT=0." >&2
    exit 1
  fi
  LOCK_FILE="${OUT_CSV}.lock"
  exec 9>"$LOCK_FILE"
  flock 9
fi

export SCRIPT_NAME INST_FOLDER INSTANCE_FROM INSTANCE_TO FAIRNESS_SET TREE_SET
export J_SET THETA_SET AVG_D_SET DEV_D_SET OUT_CSV OUT_CSV_HOUSE LOCK_OUTPUT

julia --quiet --startup-file=no --history-file=no - <<'JULIA'
include("multi.jl");
using Dates;

split_nonempty(s, sep) = [x for x in split(s, sep) if !isempty(strip(x))];
parse_int_list(s) = [parse(Int, strip(x)) for x in split_nonempty(s, ",")];
parse_float_list(s) = [parse(Float64, strip(x)) for x in split_nonempty(s, ",")];

script_name = ENV["SCRIPT_NAME"];
inst_folder = ENV["INST_FOLDER"];
instance_from = parse(Int, ENV["INSTANCE_FROM"]);
instance_to = parse(Int, ENV["INSTANCE_TO"]);
fairness_set = [String(strip(x)) for x in split_nonempty(ENV["FAIRNESS_SET"], ",")];
tree_specs = split_nonempty(ENV["TREE_SET"], ";");
J_set = parse_int_list(ENV["J_SET"]);
theta_set = parse_float_list(ENV["THETA_SET"]);
avg_d_set = parse_float_list(ENV["AVG_D_SET"]);
dev_d_set = parse_float_list(ENV["DEV_D_SET"]);
out_csv = ENV["OUT_CSV"];
out_csv_house = ENV["OUT_CSV_HOUSE"];
run_id = Dates.format(now(), "yyyymmdd_HHMMSS");

isdir(inst_folder) || error("No existe carpeta de instancias: $inst_folder");

tree_set = NamedTuple[];
for spec in tree_specs
    parts = split(spec, ":");
    length(parts) == 3 || error("TREE_SET inválido en '$spec'. Usa NBstage:childs:periods");
    push!(tree_set, (
        NBstage = parse(Int, strip(parts[1])),
        childs = parse(Int, strip(parts[2])),
        periods = parse(Int, strip(parts[3]))
    ));
end

function file_index(name::String)
    m = match(r"(\d+)", name)
    return m === nothing ? typemax(Int) : parse(Int, m.captures[1])
end
files = sort(readdir(inst_folder), by = f -> (replace(f, r"\d+" => ""), file_index(f), f));
from = max(1, instance_from);
to = min(length(files), instance_to);
from <= to || error("Rango inválido: [$instance_from,$instance_to] para $(length(files)) archivos");

configs = NamedTuple[];
for (local_file_idx, file) in enumerate(files[from:to])
    file_idx = from + local_file_idx - 1
    inFile = joinpath(inst_folder, file);
    for fair in fairness_set
        for tree in tree_set
            for J in J_set
                for theta in theta_set
                    for avg_d in avg_d_set
                        for dev_d in dev_d_set
                            push!(configs, (
                                run_id = run_id,
                                script_name = script_name,
                                inst_folder = inst_folder,
                                instance_from = instance_from,
                                instance_to = instance_to,
                                instance_no = file_idx,
                                instance_file = file,
                                instance_path = inFile,
                                fairness = fair,
                                inFile = inFile,
                                NBstage = tree.NBstage,
                                childs = tree.childs,
                                periods = tree.periods,
                                J = J,
                                theta = theta,
                                avg_d = avg_d,
                                dev_d = dev_d
                            ));
                        end
                    end
                end
            end
        end
    end
end

println("Total configuraciones a ejecutar: ", length(configs));

function load_existing_df(path::String)
    if isfile(path) && filesize(path) > 0
        try
            return CSV.read(path, DataFrame)
        catch err
            println("Warning: no se pudo leer '$path' ($(err)). Se creará un DataFrame nuevo.")
            return DataFrame()
        end
    end
    return DataFrame()
end

function reorder_columns!(df::DataFrame, front_cols::Vector{String})
    if nrow(df) == 0 && ncol(df) == 0
        return df
    end
    present_front = [c for c in front_cols if c in names(df)]
    tail = [c for c in names(df) if !(c in present_front)]
    select!(df, vcat(present_front, tail))
    return df
end

function drop_columns!(df::DataFrame, cols_to_drop::Vector{String})
    if ncol(df) == 0
        return df
    end
    keep = [c for c in names(df) if !(c in cols_to_drop)]
    select!(df, keep)
    return df
end

all_df = load_existing_df(out_csv);
all_df_house = load_existing_df(out_csv_house);
drop_cols = ["InstanceFrom", "InstanceTo", "InstanceID", "ScriptName", "RunID", "ConfigID", "InstancePath"];
drop_columns!(all_df, drop_cols);
drop_columns!(all_df_house, drop_cols);

n_total = length(configs);

summary_front = [
    "InstanceNo", "InstanceFile",
    "NBstage", "Childs", "Periods", "Theta", "Avg_d", "Dev_d", "J",
    "Fairness",
    "InstFolder"
]
house_front = [
    "InstanceNo", "InstanceFile",
    "NBstage", "Childs", "Periods", "Theta", "Avg_d", "Dev_d", "J",
    "Fairness", "House",
    "InstFolder"
]

for (cfg_idx, cfg) in enumerate(configs)
    println("Running config ", cfg_idx, "/", n_total, " -> ",
        "(fairness=", cfg.fairness,
        ", file=", basename(cfg.inFile),
        ", S=", cfg.NBstage,
        ", C=", cfg.childs,
        ", P=", cfg.periods,
        ", J=", cfg.J,
        ", theta=", cfg.theta, ")");

    df, df_house = run_single_fairness_instance(
        fairness=cfg.fairness,
        inFile=cfg.inFile,
        NBstage=cfg.NBstage,
        childs=cfg.childs,
        periods=cfg.periods,
        J=cfg.J,
        theta=cfg.theta,
        avg_d=cfg.avg_d,
        dev_d=cfg.dev_d,
        write_csv=false
    );

    df[!,"InstFolder"] = fill(cfg.inst_folder, nrow(df));
    df[!,"InstanceNo"] = fill(cfg.instance_no, nrow(df));
    df[!,"InstanceFile"] = fill(cfg.instance_file, nrow(df));
    df[!,"Fairness"] = fill(cfg.fairness, nrow(df));
    df[!,"NBstage"] = fill(cfg.NBstage, nrow(df));
    df[!,"Childs"] = fill(cfg.childs, nrow(df));
    df[!,"Periods"] = fill(cfg.periods, nrow(df));
    df[!,"J"] = fill(cfg.J, nrow(df));
    df[!,"Theta"] = fill(cfg.theta, nrow(df));
    df[!,"Avg_d"] = fill(cfg.avg_d, nrow(df));
    df[!,"Dev_d"] = fill(cfg.dev_d, nrow(df));

    df_house[!,"InstFolder"] = fill(cfg.inst_folder, nrow(df_house));
    df_house[!,"InstanceNo"] = fill(cfg.instance_no, nrow(df_house));
    df_house[!,"InstanceFile"] = fill(cfg.instance_file, nrow(df_house));
    df_house[!,"Fairness"] = fill(cfg.fairness, nrow(df_house));
    df_house[!,"NBstage"] = fill(cfg.NBstage, nrow(df_house));
    df_house[!,"Childs"] = fill(cfg.childs, nrow(df_house));
    df_house[!,"Periods"] = fill(cfg.periods, nrow(df_house));
    df_house[!,"J"] = fill(cfg.J, nrow(df_house));
    df_house[!,"Theta"] = fill(cfg.theta, nrow(df_house));
    df_house[!,"Avg_d"] = fill(cfg.avg_d, nrow(df_house));
    df_house[!,"Dev_d"] = fill(cfg.dev_d, nrow(df_house));

    drop_columns!(df, drop_cols);
    drop_columns!(df_house, drop_cols);
    reorder_columns!(df, summary_front);
    reorder_columns!(df_house, house_front);

    append!(all_df, df, cols=:union);
    append!(all_df_house, df_house, cols=:union);
    drop_columns!(all_df, drop_cols);
    drop_columns!(all_df_house, drop_cols);
    reorder_columns!(all_df, summary_front);
    reorder_columns!(all_df_house, house_front);

    # Persistir resultados en cada corrida individual.
    CSV.write(out_csv, all_df);
    CSV.write(out_csv_house, all_df_house);
end

println("Reporte resumen: ", out_csv);
println("Reporte por casa: ", out_csv_house);
JULIA

import argparse
from subprocess import run


def run_haplarithmisis(run_simg, goal, script_path, config_file, add_args=None, err_file="error.err"):
    
    r_script = f"{script_path}/{goal}/{goal}.R"
    err_file = f"{goal}_{err_file}"

    add_args_str = " ".join(add_args) if add_args else ""
    cmd = f"{run_simg} {r_script} {config_file} {add_args_str} {err_file}"
        
    result = run(cmd, capture_output=True, text=True, shell=True)
    
    print(f"Output for {goal}:\n{result.stdout}")
    if result.stderr:
        print(f"Error for {goal}:\n{result.stderr}")
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Haplarithmisis pipeline steps.")
    parser.add_argument("--run_simg", help="Singularity image command with Rscript")
    parser.add_argument("--goal", help="Pipeline step to execute")
    parser.add_argument("--script_path", help="R scripts folder")
    parser.add_argument("--config_file", help="Configuration file path")
    parser.add_argument("--add_args", nargs="*", help="Additional arguments to pass to the R script", default=[])
    parser.add_argument("--err_file", help="Output error file", default="error.err")
    args = parser.parse_args()

    run_haplarithmisis(args.run_simg, args.goal, args.script_path, args.config_file, args.add_args, args.err_file)
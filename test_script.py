import csv
import subprocess


def run_lcs_commands(input_csv, output_file):

    with open(output_file, mode="w") as out_file:

        with open(input_csv, mode="r", newline="") as file:
            reader = csv.reader(file)

            for idx, row in enumerate(reader, start=1):
                if len(row) < 2:
                    out_file.write(f"Skipping invalid row {idx}: {row}\n")
                    continue

                x = row[0]
                y = row[1]

                out_file.write(f"Running commands for pair {idx}: x = {x}, y = {y}\n")

                # serial
                run_command(["./lcs_serial", "-t", "1", "-x", x, "-y", y], out_file)

                # parallel with 4 threads
                run_command(["./lcs_parallel", "-t", "4", "-x", x, "-y", y], out_file)

                # distributed
                # runcommand(["./lcs_distributed", "-x", x, "-y", y], out_file)


def run_command(command, out_file):

    try:
        out_file.write(f"Running command: {' '.join(command)}\n")
        result = subprocess.run(command, capture_output=True, text=True, check=True)

        # Save standard output to the file
        out_file.write(f"Output:\n{result.stdout}\n")
        if result.stderr:
            out_file.write(f"Error:\n{result.stderr}\n")
    except subprocess.CalledProcessError as e:
        out_file.write(f"Error running command {' '.join(command)}: {e}\n")
        out_file.write(f"Return Code: {e.returncode}\n")
        out_file.write(f"Error Output:\n{e.stderr}\n")


input_csv = "input.csv"
output_file = "lcs_output.txt"
run_lcs_commands(input_csv, output_file)

print(f"Specific random pairs and their outputs have been saved to '{output_file}'.")

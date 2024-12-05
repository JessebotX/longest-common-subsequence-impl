import random
import string
import csv


# Generate a random string of the given length.
def generate_random_string(length):
    return "".join(random.choices(string.ascii_letters, k=length))


# Generate pairs of random strings with specific lengths.
def generate_specific_pairs():

    pairs = [
        (1, 1),
        (10000, 1),
        (1, 100000),
        (10, 10),
        (100, 100),
        (1000, 1000),
        (10000, 10000),
        (100000, 100000),
    ]

    random_pairs = []
    for x, y in pairs:
        str1 = generate_random_string(x)
        str2 = generate_random_string(y)
        random_pairs.append((str1, str2))

    return random_pairs


# Save the pairs of strings to a CSV file.
def save_pairs_to_csv(pairs, filename="input.csv"):

    with open(filename, mode="w", newline="") as file:
        writer = csv.writer(file)
        for str1, str2 in pairs:
            writer.writerow([str1, str2])


# Generate specific pairs
pairs = generate_specific_pairs()

# Save to CSV file
save_pairs_to_csv(pairs)

print("Specific random pairs have been saved to 'input.csv'.")

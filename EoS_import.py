def file_read(input_file):
    data_list = []
    with open(input_file, newline='\n') as csvfile:
        file_read = csv.reader(csvfile, delimiter=',', quotechar='|')
        data_list = [row for row in file_read]
    pressure_list = [el for el in data_list[el][0]]
    radius_list = [el for el in data_list[el][1]]

    return pressure_list, radius_list
def EOS_import(input_file):
    pressure, radius = file_read(input_file)

    return pressure,radius
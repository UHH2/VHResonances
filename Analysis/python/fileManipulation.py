import sys
import os
import subprocess
import time

def is_commented(line, filename):
    #ignore empty line
    if len(line.split())<1: return False
    #get rid of whitespace
    beginning = line.split()[0]
    if ".c" in filename:
        if beginning[0:2] == "//": return True
        else: return False
    elif ".py" in filename:
        if beginning[0:1] == "#": return True
        else: return False
    elif ".sh" in filename:
        if beginning[0:1] == "#": return True
        else: return False
    elif ".xml" in filename:
        if beginning[0:4] == "<!--": return True
        else: return False
    else:
        print "File extension not known!"
        return False

def number_leading_whitespace(line):
    n = 0
    while line[n] == " ":
        n += 1
    return n

def get_whitespace(n):
    line = ""
    for i in range(0, n):
        line += " "
    return line

# return the commented line
def comment_single_line(filename, line):
    if is_commented(line, filename):
        print "Warning: This line is already commented."
        return line
    elif len(line.split())<1:
        print "Warning: This line is empty."
        return line
    else:
        n = number_leading_whitespace(line)
        startline = ""
        for i in range(0, n):
            startline += " "
        if ".c" in filename: return     startline + "// "   + line
        elif ".py" in filename: return  startline + "# "    + line
        elif ".sh" in filename: return  startline + "# "    + line
        elif ".xml" in filename: return startline + "<!-- " + line[n:len(line)-1] + " -->\n"

# return the uncommented line
def uncomment_single_line(filename, line):
    if not is_commented(line, filename):
        print "Warning: This line is already uncommented."
        return line
    else:
        if ".c" in filename:
            n = number_leading_whitespace(line)
            newline = line[0:n]
            newline += line[n+3:]
            return newline
        elif ".py" in filename:
            n = number_leading_whitespace(line)
            newline = line[0:n]
            newline += line[n+2:]
            return newline
        elif ".xml" in filename:
            n = number_leading_whitespace(line)
            newline = line[0:n]
            newline += line[n+5:len(line)-5]
            return newline + "\n"

# comment all lines that match the single control
def comment_line(path, filename, controls, comment=True):
    newText = []
    with open(path+filename, "U") as file:
        lines = file.readlines()
        for line in lines:
            isToChange = True
            for control in controls:
                isToChange = isToChange and control in line
            if isToChange:
                if comment:
                    newLine = comment_single_line(filename, line)
                else:
                    newLine = uncomment_single_line(filename, line)
                newText.append(newLine)
            else:
                newText.append(line)
    with open(path+filename, "w") as outputfile:
        for line in newText:
            outputfile.write(line)

# comment all lines that match any of the controls
def comment_lines(path, filename, controls, comment=True, remove=False):
    newText = []
    with open(path+filename, "U") as file:
        lines = file.readlines()
        for line in lines:
            isToSave = True
            for control in controls:
                isToChange = True
                for ElToCheck in control:
                    isToChange = isToChange and ElToCheck in line
                if isToChange:
                    isToSave = False
                    if comment:
                        newLine = comment_single_line(filename, line)
                    else:
                        newLine = uncomment_single_line(filename, line)
                    if not remove:
                        newText.append(newLine)
            if isToSave:
                newText.append(line)
    with open(path+filename, "w") as outputfile:
        for line in newText:
            outputfile.write(line)


# replace input to outputs in the line that match any of the controls
def change_lines(path, filename, controls, inputs, outputs):
    newText = []
    with open(path+filename, "U") as file:
        lines = file.readlines()
        for line in lines:
            isToSave = True
            for index, control in enumerate(controls):
                isToChange = True
                for ElToCheck in control:
                    isToChange = isToChange and ElToCheck in line
                for index_1, el in enumerate(inputs[index]):
                    isToChange = isToChange and inputs[index][index_1] in line
                if isToChange:
                    isToSave = False
                    newLine = line
                    for index_1, el in enumerate(inputs[index]):
                        newLine = newLine.replace(inputs[index][index_1],outputs[index][index_1])
                    newText.append(newLine)
            if isToSave:
                newText.append(line)
    with open(path+filename, "w") as outputfile:
        for line in newText:
            outputfile.write(line)


# replace only one line
def change_line(path, filename, controls, inputs, outputs):
    newText = []
    with open(path+filename, "U") as file:
        lines = file.readlines()
        for line in lines:
            isToChange = True
            for control in controls:
                isToChange = isToChange and control in line
            if isToChange:
                newLine = line
                for index, el in enumerate(inputs):
                    newLine = newLine.replace(inputs[index],outputs[index])
                newText.append(newLine)
            else:
                newText.append(line)
    with open(path+filename, "w") as outputfile:
        for line in newText:
            outputfile.write(line)

# it's slower wrt change_line
def sub_line(path, filename, controls, inputs, outputs):
    newText = []
    with open(path+filename, "U") as file:
        lines = file.readlines()
        for i, line in enumerate(lines):
            isToChange = True
            for control in controls:
                isToChange = isToChange and control in line
            if isToChange:
                for index, el in enumerate(inputs):
                    command = ["sed", "-i -e"]
                    command.append(str(index)+"s/"+inputs[index]+"/"+outputs[index]+"/g")
                    command.append(path+filename)
                    process = subprocess.Popen(command)
                    process.wait()

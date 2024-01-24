from tkinter import * # used for rendering graphics on a display window
from gtts import gTTS # used for converting text to voice -google text to speech-
from playsound import playsound # used for playing the converter voice from the text

        ### Initialized window ###
# Tk() to initialized tkinter which will be used for GUI
# geometry() used to set the width and height of the window
# configure() used to access window attributes
# bg will used to set the color of the background
# title() set the title of the window

root = Tk()
root.geometry('500x250')
root.resizable(0,0)
root.config(bg = 'purple')
root.title('TEXT2SPEECH')

# Label() widget is used to display one or more than one line of text that users can’t able to modify.

# root is the name which we refer to our window
# text which we display on the label
# font in which the text is written
# pack organized widget in block
# Msg is a string type variable
# Entry() used to create an input text field
# textvariable used to retrieve the current text to entry widget
# place() organizes widgets by placing them in a specific position in the parent widget

#heading
Label(root, text = 'TEXT2SPEECH' , font='arial 30 bold' , bg ='orange').pack()

#label
Label(root, text ='Enter Text', font ='arial 15 bold', bg ='orange').place(x=20,y=60)

#text variable
Msg = StringVar()

#Entry
entry_field = Entry(root,textvariable =Msg, width ='50')
entry_field.place(x=20 , y=100)


        ### Function to Convert Text to Speech ###

def Text_to_speech():
    Message = entry_field.get()
    speech = gTTS(text = Message)
    speech.save('audio_file.mp3')
    playsound('audio_file.mp3')

# Message variable will stores the value of entry_field
# text is the sentences or text to be read.
# lang takes the language to read the text. The default language is English.
# slow use to reads text more slowly. The default is False.
# As we want the default value of lang, so no need to give that to gTTS.

# speech stores the converted voice from the text
# speech.save will saves the converted file  as mp3 file
# playsound() used to play the sound

def Exit():
    root.destroy()
# root.destroy() will quit the program by stopping the mainloop()

def Reset():
    Msg.set("")
# reset function set Msg variable to empty strings

#Button
Button(root, text = "PLAY" , font = 'arial 15 bold', command = Text_to_speech, width =4).place(x=25, y=140)
Button(root,text = 'EXIT',font = 'arial 15 bold' , command = Exit, bg = 'OrangeRed1').place(x=100,y=140)
Button(root, text = 'RESET', font='arial 15 bold', command = Reset).place(x=175 , y =140)
# Button() widget used to display button on the window

#infinite loop to run program
root.mainloop()
# root.mainloop() is a method that executes when we want to run our program

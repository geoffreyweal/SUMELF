"""
folder_methods.py, Geoffrey Weal, 19/7/22

This script is designed to hold methods useful for creating and removing folders.
"""
import os, shutil

def make_folder(folder_path):
	"""
	This method will create the folder given by folder_path if it does not exist.

	Parameters
	----------
	folder_path : str. 
		The path to the folder you want to create	
	"""
	if not os.path.exists(folder_path):
		os.makedirs(folder_path)

def remove_folder(folder_path):
	"""
	This method will remove the folder given by folder_path

	Parameters
	----------
	folder_path : str. 
		The path to the folder you want to remove.	
	"""
	if os.path.exists(folder_path):
		shutil.rmtree(folder_path)

def move_folder(from_folder_path, to_folder_path):
	"""
	This method will move a folder from from_folder_path to to_folder_path

	Parameters
	----------
	from_folder_path : str. 
		The path of the folder to move
	from_folder_path : str. 
		The path to move the from_folder_path folder into 
	"""
	shutil.move(from_folder_path,to_folder_path)

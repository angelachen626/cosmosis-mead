import ctypes as ct

dll = ct.cdll.LoadLibrary("libcosmosis.so")

c_block = ct.c_size_t
c_status = ct.c_int
c_str = ct.c_char_p
c_int = ct.c_int
c_int_p = ct.POINTER(ct.c_int)

def load_library_function(namespace, name, argtypes, restype):
	function = getattr(dll,name)
	function.argtypes = argtypes
	function.restype = restype
	namespace[name] = function

class c_complex(ct.Structure):
	_fields_ = [("real", ct.c_double), ("imag", ct.c_double)]


def load_function_types(namespace, c_type, c_name):
	load_library_function(namespace, "c_datablock_put_%s"%c_name, [c_block, c_str, c_str, c_type], c_status)
	load_library_function(namespace, "c_datablock_replace_%s"%c_name, [c_block, c_str, c_str, c_type], c_status)
	load_library_function(namespace, "c_datablock_get_%s"%c_name, [c_block, c_str, c_str, ct.POINTER(c_type)], c_status)

def load_array_function_types(namespace, c_type, c_name):
	load_library_function(namespace, "c_datablock_put_%s_array_1d"%c_name, [c_block, c_str, c_str, ct.POINTER(c_type), c_int], c_status)
	load_library_function(namespace, "c_datablock_get_%s_array_1d"%c_name, [c_block, c_str, c_str, ct.POINTER(ct.POINTER(c_type)), c_int_p], c_status)
	load_library_function(namespace, "c_datablock_get_%s_array_1d_preallocated"%c_name, [c_block, c_str, c_str, ct.POINTER(c_type), c_int_p, c_int], c_status)



load_function_types(locals(), ct.c_int, 'int')
load_function_types(locals(), ct.c_double, 'double')
load_function_types(locals(), c_complex, 'complex')

load_array_function_types(locals(), ct.c_int, 'int')

load_library_function(
	locals(), 
	"make_c_datablock",
	[],
	c_block
)

load_library_function(
	locals(), 
	"make_c_datablock",
	[],
	c_block
)




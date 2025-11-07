import ROOT
import sys

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <ntuple_name> <file.root>")
        sys.exit(1)

    ntuple_name = sys.argv[1]
    file_name   = sys.argv[2]

    reader = ROOT.RNTupleReader.Open(ntuple_name, file_name)

    # Print summary info: schema, number of entries, compression, etc.
    reader.PrintInfo(ROOT.ENTupleInfo.kStorageDetails)
    reader.PrintInfo(ROOT.ENTupleInfo.kSummary)

if __name__ == "__main__":
    main()


# #!/usr/bin/env python3
# """
# Script to dump all RNTuple fields from a given collection in a ROOT file using PyROOT.
# """

# import ROOT
# import argparse
# import sys
# from typing import Optional


# def dump_rntuple_fields(file_path: str, collection_name: str, max_entries: Optional[int] = None) -> None:
#     """
#     Dump all fields from an RNTuple collection in a ROOT file using PyROOT.
    
#     Args:
#         file_path: Path to the ROOT file
#         collection_name: Name of the RNTuple collection
#         max_entries: Maximum number of entries to display (None for all)
#     """
#     # try:
#     #     # Open the ROOT file
#     #     root_file = ROOT.TFile.Open(file_path, "READ")
#     #     if not root_file or root_file.IsZombie():
#     #         print(f"Error: Cannot open file '{file_path}'")
#     #         return
        
#     #     print(f"Opened ROOT file: {file_path}")
        
#     #     # List all keys in the file
#     #     keys = []
#     #     key_list = root_file.GetListOfKeys()
#     #     for i in range(key_list.GetSize()):
#     #         keys.append(key_list.At(i).GetName())
#     #     print(f"Available keys: {keys}")
        
#     # Try to open the RNTuple
#     try:
#         # Create RNTupleReader
#         reader = ROOT.RNTupleReader.Open(collection_name, file_path)
#         print(f"\nSuccessfully opened RNTuple: {collection_name}")
        
#     except Exception as e:
#         print(f"Error: Cannot open RNTuple '{collection_name}': {e}")
#         print(f"Available keys: {keys}")
#         root_file.Close()
#         return
    
#     # Get descriptor to access field information
#     descriptor = reader.GetDescriptor()
    
#     # Get total number of entries
#     total_entries = reader.GetNEntries()
#     print(f"Total entries: {total_entries}")
    
#     # Get field information
#     field_range = descriptor.GetFieldRange()
#     fields_info = []
    
#     print(f"\nField Information:")
#     print("=" * 80)
    
#     # Iterate through fields
#     for field_id in field_range:
#         field = descriptor.GetFieldDescriptor(field_id)
#         field_name = field.GetFieldName()
#         field_type = field.GetTypeName()
#         parent_id = field.GetParentId()
        
#         # Only show top-level fields (no parent or parent is root)
#         if parent_id == ROOT.Experimental.kInvalidDescriptorId or parent_id == 0:
#             fields_info.append((field_name, field_type, field_id))
#             print(f"Field: {field_name}")
#             print(f"  Type: {field_type}")
#             print(f"  Field ID: {field_id}")
#             print()
    
#     # Determine how many entries to show
#     entries_to_show = total_entries
#     if max_entries is not None:
#         entries_to_show = min(max_entries, total_entries)
    
#     if entries_to_show > 1000:
#         entries_to_show = 1000
#         print(f"Limiting display to first {entries_to_show} entries for performance")
    
#     print(f"Showing first {entries_to_show} entries:")
#     print("=" * 120)
    
#     # Create value objects for each field
#     field_values = {}
#     for field_name, field_type, field_id in fields_info:
#         try:
#             # Create appropriate value object based on type
#             if "int" in field_type.lower() or field_type in ["Int_t", "UInt_t"]:
#                 field_values[field_name] = reader.GetModel().MakeField["int"](field_name)
#             elif "float" in field_type.lower() or field_type in ["Float_t", "Double_t"]:
#                 field_values[field_name] = reader.GetModel().MakeField["float"](field_name)
#             elif "double" in field_type.lower():
#                 field_values[field_name] = reader.GetModel().MakeField["double"](field_name)
#             elif "bool" in field_type.lower():
#                 field_values[field_name] = reader.GetModel().MakeField["bool"](field_name)
#             elif "string" in field_type.lower() or field_type == "std::string":
#                 field_values[field_name] = reader.GetModel().MakeField["std::string"](field_name)
#             else:
#                 # Try generic approach
#                 try:
#                     field_values[field_name] = reader.GetModel().MakeField[field_type](field_name)
#                 except:
#                     print(f"Warning: Cannot create value object for field '{field_name}' of type '{field_type}'")
#                     continue
                    
#             # Connect to reader
#             reader.GetModel().ConnectModel(field_values[field_name])
            
#         except Exception as e:
#             print(f"Warning: Error setting up field '{field_name}': {e}")
#             continue
    
#     # Print header
#     header = "Entry".ljust(8)
#     valid_fields = []
#     for field_name, _, _ in fields_info:
#         if field_name in field_values:
#             header += f"{field_name}"[:20].ljust(22)
#             valid_fields.append(field_name)
#     print(header)
#     print("-" * len(header))
    
#     # Read and display entries
#     display_entries = min(20, entries_to_show)  # Limit to 20 for readability
    
#     for entry_idx in range(display_entries):
#         try:
#             reader.LoadEntry(entry_idx)
            
#             row = f"{entry_idx}".ljust(8)
#             for field_name in valid_fields:
#                 try:
#                     value = str(field_values[field_name].GetValue())[:20]
#                     row += f"{value}".ljust(22)
#                 except Exception as e:
#                     row += "ERROR".ljust(22)
            
#             print(row)
            
#         except Exception as e:
#             print(f"Error reading entry {entry_idx}: {e}")
#             break
    
#     if display_entries < entries_to_show:
#         print(f"\n... (showing only first {display_entries} entries for readability)")
    
#     # Show some statistics for numeric fields
#     print(f"\nField Statistics (first {entries_to_show} entries):")
#     print("=" * 80)
    
#     for field_name in valid_fields:
#         try:
#             values = []
#             for entry_idx in range(min(entries_to_show, 1000)):  # Limit for performance
#                 reader.LoadEntry(entry_idx)
#                 try:
#                     val = field_values[field_name].GetValue()
#                     if isinstance(val, (int, float)):
#                         values.append(val)
#                 except:
#                     continue
            
#             if values:
#                 print(f"{field_name}:")
#                 print(f"  Count: {len(values)}")
#                 print(f"  Min: {min(values)}")
#                 print(f"  Max: {max(values)}")
#                 print(f"  Mean: {sum(values)/len(values):.3f}")
#                 print()
                
#         except Exception as e:
#             continue
        
#     # Clean up
#     root_file.Close()
        
#     # except Exception as e:
#     #     print(f"Error: {e}")
#     #     import traceback
#     #     traceback.print_exc()


# def list_collections(file_path: str) -> None:
#     """List all available collections in the ROOT file."""
#     try:
#         root_file = ROOT.TFile.Open(file_path, "READ")
#         if not root_file or root_file.IsZombie():
#             print(f"Error: Cannot open file '{file_path}'")
#             return
        
#         print(f"Objects in {file_path}:")
        
#         # List all keys
#         key_list = root_file.GetListOfKeys()
#         for i in range(key_list.GetSize()):
#             key = key_list.At(i)
#             obj_name = key.GetName()
#             obj_class = key.GetClassName()
            
#             # Try to determine if it's an RNTuple
#             is_rntuple = False
#             try:
#                 # Try to open as RNTuple
#                 test_reader = ROOT.Experimental.RNTupleReader.Open(obj_name, file_path)
#                 is_rntuple = True
#                 test_reader = None  # Clean up
#             except:
#                 pass
            
#             rntuple_indicator = " (RNTuple)" if is_rntuple else ""
#             print(f"  {obj_name} ({obj_class}){rntuple_indicator}")
        
#         root_file.Close()
        
#     except Exception as e:
#         print(f"Error listing collections: {e}")


# def main():
#     parser = argparse.ArgumentParser(
#         description="Dump RNTuple fields from a ROOT file using PyROOT",
#         formatter_class=argparse.RawDescriptionHelpFormatter,
#         epilog="""
# Examples:
#   %(prog)s myfile.root MyRNTuple
#   %(prog)s myfile.root MyRNTuple --max-entries 100
#   %(prog)s --list myfile.root

# Note: Requires ROOT with RNTuple support and PyROOT installed.
#         """
#     )
    
#     parser.add_argument("file", help="Path to the ROOT file")
#     parser.add_argument("collection", nargs='?', help="Name of the RNTuple collection")
#     parser.add_argument("--max-entries", "-n", type=int, 
#                        help="Maximum number of entries to display")
#     parser.add_argument("--list", "-l", action="store_true",
#                        help="List all collections in the file")
    
#     args = parser.parse_args()
    
#     # Enable RNTuple experimental features
#     try:
#         ROOT.ROOT.EnableImplicitMT()
#         ROOT.gSystem.Load("libROOTNTuple")
#     except:
#         print("Warning: Could not load RNTuple libraries. Make sure ROOT is built with RNTuple support.")
    
#     if args.list:
#         list_collections(args.file)
#     elif args.collection:
#         dump_rntuple_fields(args.file, args.collection, args.max_entries)
#     else:
#         print("Error: Please specify a collection name or use --list to see available collections.")
#         parser.print_help()
#         sys.exit(1)


# if __name__ == "__main__":
#     main()
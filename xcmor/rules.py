from .log import get_logger
from .mapping import dtype_map

logger = get_logger(__name__)


class rules:
    drop = True

    @classmethod
    def type(cls, obj):
        """apply dtype rule to dataarray"""
        dtype = obj.attrs["type"]
        if obj.dtype != dtype_map[dtype]:
            logger.warning(
                #   f"converting {obj.name or "data"} from "
                f"{obj.dtype} to {dtype_map[dtype]}"
            )
            obj = obj.astype(dtype_map[dtype])
        if cls.drop is True:
            del obj.attrs["type"]
        return obj

    @classmethod
    def out_name(cls, obj):
        """apply renaming rule"""
        out_name = obj.attrs["out_name"]
        obj.name = out_name
        # obj.rename_vars({v: out_name})
        if cls.drop is True:
            del obj.attrs["out_name"]
        return obj

    @classmethod
    def valid_min(cls, obj):
        valid_min = obj.attrs.get("valid_min")
        if valid_min:
            try:
                assert obj.min() >= valid_min
            except Exception:
                raise Exception(
                    f"{obj.name or 'data'} is violating valid_min: {valid_min}"
                )
        if cls.drop is True:
            del obj.attrs["valid_min"]
        return obj

    @classmethod
    def valid_max(cls, obj):
        valid_max = obj.attrs.get("valid_max")
        if valid_max:
            try:
                assert obj.min() <= valid_max
            except Exception:
                raise Exception(
                    f"{obj.name or 'data'} is violating valid_min: {valid_max}"
                )
        if cls.drop is True:
            del obj.attrs["valid_max"]
        return obj

    @classmethod
    def standard_name(cls, obj):
        from cf_xarray.utils import parse_cf_standard_name_table

        sname = obj.attrs["standard_name"]
        info, table, aliases = parse_cf_standard_name_table()
        if sname not in table.keys():
            raise Exception(f"{sname} is not a valid standard name.")
        return obj
